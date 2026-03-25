function [S, state, normalized_output_frame, tau_map] = adaptiveTimeSurfaceZhu(...
    x_list, y_list, t_list, p_list, imgSz, state, params)
% ADAPTIVETIMESURFACEZHU  Adaptive Time Surface (Zhu et al., IROS 2023).
%
%   [S, STATE, NORMALIZED_OUTPUT_FRAME, TAU_MAP] = adaptiveTimeSurfaceZhu(...
%       X, Y, T, P, IMGSZ, STATE, PARAMS)
%
%   Reference implementation of the ATS algorithm from:
%       Zhu, S., Tang, Z., Yang, M., Learned-Miller, E., and Kim, D.,
%       "Event Camera-based Visual Odometry for Dynamic Motion Tracking
%       of a Legged Robot Using Adaptive Time Surface," 2023 IEEE/RSJ
%       International Conference on Intelligent Robots and Systems
%       (IROS), Detroit, MI, Oct. 1-5, 2023, pp. 3475-3482.
%       DOI: 10.1109/IROS55552.2023.10342048
%
%   Zhu's ATS computes a pixel-wise decay rate from the average temporal
%   staleness of neighboring pixels within a fixed 16-pixel spatial
%   pattern. The surface itself is the classical Lagorce et al. (2017)
%   exponential time surface with the spatially varying tau substituted
%   for the fixed constant.
%
%   Inputs:
%     X, Y    - [N x 1] vectors of pixel coordinates (row, col; 1-indexed)
%     T       - [N x 1] vector of timestamps [seconds]
%     P       - [N x 1] vector of polarities (+1 or -1)  [unused by
%               the algorithm but accepted for API consistency]
%     IMGSZ   - [1 x 2] image dimensions [nRows, nCols]
%     STATE   - Struct with persistent memory fields:
%                 .t_last  [imgSz] - timestamp of most recent event
%     PARAMS  - Struct with ATS parameters:
%                 .tau_u      - Upper bound on decay rate [s] (default: 0.5)
%                 .tau_l      - Lower bound on decay rate [s] (default: 0.01)
%                 .n_neighbors - Number of most-recent neighbors to use
%                                in the staleness average (default: 8)
%                 .blur_sigma  - Gaussian blur sigma for post-smoothing
%                                of the output surface [pixels] (default: 1.0)
%                 .median_sz   - Median filter kernel size [pixels, odd]
%                                (default: 3)
%
%   Outputs:
%     S                       - Raw ATS surface [imgSz], values in [0, 1]
%     STATE                   - Updated state struct
%     NORMALIZED_OUTPUT_FRAME - Display-ready surface mapped to [0, 1]
%                               after blur + median filter (as in paper)
%     TAU_MAP                 - [imgSz] Per-pixel adaptive decay rate
%
%   Algorithm (Eqs. 1-2 of the paper):
%     1. Update t_last map with incoming events.
%     2. For each pixel, sample 16 neighbors in a fixed spatial pattern
%        (Fig. 3 of paper). Rank by recency, select the n most recent.
%     3. Compute adaptive tau:
%          tau(x) = max(tau_u - (1/n) * sum(t - t_last_i), tau_l)
%     4. Compute time surface:
%          T(x,t) = exp(-(t - t_last(x)) / tau(x))
%     5. Apply Gaussian blur + median blur for smoothing.
%
%   Notes:
%     - This is a batch-mode reference implementation suitable for the
%       IEI-ATS frame-based processing pipeline. The original paper
%       does not describe an optimized implementation.
%     - Events must be pre-sorted in ascending temporal order.
%     - Polarity values are accepted for API consistency but are not
%       used by the algorithm (the surface is unsigned).
%     - The 16-pixel spatial pattern is reconstructed from Fig. 3 of
%       the paper. It samples pixels at distances of approximately
%       3-4 pixels in a roughly circular arrangement.
%     - SIGN CONVENTION NOTE: Eq. 2 as published produces larger tau
%       (slower decay) when neighbors are fresh (small staleness) and
%       smaller tau (faster decay) when neighbors are stale. The
%       paper's text describes the opposite behavior. This
%       implementation follows the equation as published. See the
%       companion analysis document for discussion.
%
%   Structural comparison with IEI-ATS:
%     - Reset-based surface (no recursive state; t_last overwritten)
%     - Neighborhood staleness drives tau (not per-pixel IEI)
%     - Unsigned output (no polarity preservation)
%     - No pre-filtering, no attack-release, no divisive normalization
%
%   See also: accumulator.localAdaptiveTimeSurface,
%             accumulator.motionEncodedTimeSurface,
%             accumulator.timeSurface

    % ----------------------------------------------------------------
    % 0. Parse parameters with defaults
    % ----------------------------------------------------------------
    tau_u       = params.tau_u;          % Upper bound on decay [s]
    tau_l       = params.tau_l;          % Lower bound on decay [s]
    n_neighbors = params.n_neighbors;    % # of neighbors for avg
    blur_sigma  = params.blur_sigma;     % Gaussian blur sigma [px]
    median_sz   = params.median_sz;      % Median filter size [px]

    [H, W] = deal(imgSz(1), imgSz(2));

    % ----------------------------------------------------------------
    % 1. Define the 16-pixel spatial neighborhood pattern (Fig. 3)
    % ----------------------------------------------------------------
    % Reconstructed from Fig. 3 of the paper. The pattern samples
    % 16 pixels at approximately 3-4 pixel distances in a roughly
    % circular arrangement around the target pixel. The exact offsets
    % are not specified numerically in the paper, so we use a
    % symmetric pattern that matches the figure's geometry.
    %
    % Layout (Fig. 3): pixels at cardinal + diagonal + intermediate
    % directions at distances ~3 px, giving 16 points on a rough
    % circle of radius 3.
    %
    % Pattern (dr, dc) offsets from the target pixel:
    pattern = [
        -3,  0;   % N
        -2, -2;   % NW
         0, -3;   % W
         2, -2;   % SW
         3,  0;   % S
         2,  2;   % SE
         0,  3;   % E
        -2,  2;   % NE
        -1, -3;   % WNW upper
        -3, -1;   % NNW
        -3,  1;   % NNE
        -1,  3;   % ENE
         1,  3;   % ESE
         3,  1;   % SSE
         3, -1;   % SSW
         1, -3;   % WNW lower
    ];
    n_pattern = size(pattern, 1);  % Should be 16

    % ----------------------------------------------------------------
    % 2. Update t_last map with incoming events
    % ----------------------------------------------------------------
    x_list = x_list(:);
    y_list = y_list(:);
    t_list = t_list(:);

    if isempty(t_list)
        S = zeros(imgSz);
        normalized_output_frame = zeros(imgSz);
        tau_map = ones(imgSz) * tau_u;
        return;
    end

    % Process events sequentially (latest timestamp wins per pixel)
    for k = 1:numel(t_list)
        state.t_last(x_list(k), y_list(k)) = t_list(k);
    end

    % Current evaluation time: timestamp of the last event in batch
    t_now = t_list(end);

    % ----------------------------------------------------------------
    % 3. Compute per-pixel adaptive tau (Eq. 2)
    % ----------------------------------------------------------------
    tau_map = ones(imgSz) * tau_u;   % Default to upper bound

    % Pad t_last map to handle boundary pixels.
    % Use zero-padding (equivalent to "never activated" neighbors).
    pad = 4;  % Max absolute offset in the pattern is 3, pad by 4
    t_last_padded = zeros(H + 2*pad, W + 2*pad);
    t_last_padded(pad+1:pad+H, pad+1:pad+W) = state.t_last;

    % Vectorized: collect neighbor timestamps for all pixels at once.
    % Result shape: [H, W, n_pattern]
    neighbor_timestamps = zeros(H, W, n_pattern);
    for p_idx = 1:n_pattern
        dr = pattern(p_idx, 1);
        dc = pattern(p_idx, 2);
        neighbor_timestamps(:,:,p_idx) = ...
            t_last_padded(pad+1+dr : pad+H+dr, pad+1+dc : pad+W+dc);
    end

    % Exclude neighbors that have never fired (t_last == 0).
    % Replace zeros with NaN so they sort to the end.
    neighbor_timestamps(neighbor_timestamps == 0) = NaN;

    % Sort descending along dim 3 (most recent first).
    sorted_neighbors = sort(neighbor_timestamps, 3, 'descend', ...
        'MissingPlacement', 'last');

    % Select the n_neighbors most recent.
    n_sel = min(n_neighbors, n_pattern);
    selected = sorted_neighbors(:, :, 1:n_sel);

    % Compute staleness: (t_now - t_last_i) for each selected neighbor.
    staleness = t_now - selected;

    % Count valid (non-NaN) neighbors per pixel.
    valid_mask = ~isnan(selected);
    valid_count = sum(valid_mask, 3);

    % Sum of staleness for valid neighbors only.
    staleness(isnan(staleness)) = 0;
    staleness_sum = sum(staleness, 3);

    % Mean staleness (avoid division by zero).
    mean_staleness = zeros(imgSz);
    has_valid = valid_count > 0;
    mean_staleness(has_valid) = staleness_sum(has_valid) ./ ...
        valid_count(has_valid);

    % Eq. 2: tau(x) = max(tau_u - mean_staleness, tau_l)
    tau_map = max(tau_u - mean_staleness, tau_l);

    % Pixels with NO valid neighbors keep tau = tau_u (slow decay).
    tau_map(~has_valid) = tau_u;

    % ----------------------------------------------------------------
    % 4. Compute time surface (Eq. 1)
    % ----------------------------------------------------------------
    % T(x, t) = exp(-(t - t_last(x)) / tau(x))
    dt_map = t_now - state.t_last;

    S = exp(-dt_map ./ tau_map);

    % Zero out pixels that have never fired.
    S(state.t_last == 0) = 0;

    % Clamp to [0, 1] (should already be by construction).
    S = max(0, min(1, S));

    % ----------------------------------------------------------------
    % 5. Post-process: Gaussian blur + median blur
    %    (Section II-A: "Then blur and median blur filters are applied
    %     to produce a smoother result.")
    % ----------------------------------------------------------------
    if blur_sigma > 0
        smoothed = imgaussfilt(S, blur_sigma);
    else
        smoothed = S;
    end

    if median_sz > 1
        smoothed = medfilt2(smoothed, [median_sz, median_sz]);
    end

    normalized_output_frame = smoothed;

end