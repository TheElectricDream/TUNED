function [S, state, normalized_output_frame] = motionEncodedTimeSurface(...
    x_list, y_list, t_list, p_list, imgSz, state, params)
% MOTIONENCODEDTIMESURFACE  Motion-Encoded Time-Surface (METS).
%
%   [S, STATE, NORMALIZED_OUTPUT_FRAME] = motionEncodedTimeSurface(...
%       X, Y, T, P, IMGSZ, STATE, PARAMS)
%
%   Reference implementation of the METS algorithm from:
%       Xu et al., "METS: Motion-Encoded Time-Surface for Event-Based
%       High-Speed Pose Tracking," International Journal of Computer
%       Vision, vol. 133, pp. 4401-4419, 2025.
%       DOI: 10.1007/s11263-025-02379-6
%
%   METS dynamically encodes a pixel-wise decay rate based on the
%   estimated instantaneous edge velocity at each pixel. The velocity
%   is derived from the spatio-temporal ordering of recent event
%   timestamps within a local observation window. The decay is then
%   normalized to edge displacement rather than elapsed time, making
%   the surface invariant to motion speed.
%
%   Inputs:
%     X, Y    - [N x 1] vectors of pixel coordinates (row, col; 1-indexed)
%     T       - [N x 1] vector of timestamps [seconds]
%     P       - [N x 1] vector of polarities (+1 or -1)
%     IMGSZ   - [1 x 2] image dimensions [nRows, nCols]
%     STATE   - Struct with persistent memory fields:
%                 .t_last_pos  [imgSz] - last event timestamp (positive pol.)
%                 .t_last_neg  [imgSz] - last event timestamp (negative pol.)
%                 .t_last_any  [imgSz] - last event timestamp (either pol.)
%                 .p_last      [imgSz] - polarity of last event at pixel
%     PARAMS  - Struct with METS parameters:
%                 .R     - Observation window half-size [pixels] (default: 4)
%                 .n     - Velocity estimation range [pixels]   (default: 3)
%                 .d     - Decay step [pixels]                  (default: 5)
%                 .d_th  - Decay distance threshold [pixels]    (default: 8)
%
%   Outputs:
%     S                       - Raw METS surface [imgSz], values in [0, 1]
%     STATE                   - Updated state struct
%     NORMALIZED_OUTPUT_FRAME - Display-ready surface mapped to [0, 1]
%
%   Algorithm (Eqs. 4-7 of the paper):
%     1. For each event, update the polarity-specific timestamp map.
%     2. For each pixel that received an event in this batch, estimate
%        the local edge velocity v(x,y) from the temporal ordering of
%        events within a (2R+1)x(2R+1) observation window (Eq. 6).
%     3. Compute the METS surface at t_now using the motion-encoded
%        exponential decay (Eq. 7):
%          METS(x,y) = exp(-(t_now - t_last) / (d * v))  if active
%                    = 0                                   if beyond d_th
%
%   Notes:
%     - This is a batch-mode reference implementation suitable for the
%       IEI-ATS frame-based processing pipeline. The original paper
%       describes an optimized event-by-event C++ implementation using
%       doubly-linked lists and a look-up architecture for real-time
%       performance.
%     - Events must be pre-sorted in ascending temporal order.
%     - Polarity values must be +1 or -1 (not 0/1).

    % ----------------------------------------------------------------
    % 0. Parse parameters with defaults from Section 3.3 of the paper
    % ----------------------------------------------------------------
    R    = params.R;       % Observation window half-size (default: 4)
    n    = params.n;       % Velocity estimation range   (default: 3)
    d    = params.d;       % Decay step                  (default: 5)
    d_th = params.d_th;    % Decay distance threshold    (default: 8)

    [H, W] = deal(imgSz(1), imgSz(2));
    win_side = 2*R + 1;    % Observation window side length

    % Number of timestamps in the temporal sequence (Eq. 5)
    seq_len = n * win_side + 1;

    % ----------------------------------------------------------------
    % 1. Update timestamp maps event-by-event
    % ----------------------------------------------------------------
    % Ensure column vectors
    x_list = x_list(:);
    y_list = y_list(:);
    t_list = t_list(:);
    p_list = p_list(:);

    if isempty(t_list)
        S = zeros(imgSz);
        normalized_output_frame = zeros(imgSz);
        return;
    end

    % Track which pixels received events in this batch
    active_mask = false(imgSz);

    % Process events sequentially to maintain correct temporal ordering
    % (latest timestamp wins when multiple events hit the same pixel)
    for k = 1:numel(t_list)
        r = x_list(k);
        c = y_list(k);
        
        if p_list(k) > 0
            state.t_last_pos(r, c) = t_list(k);
        else
            state.t_last_neg(r, c) = t_list(k);
        end
        
        state.t_last_any(r, c) = t_list(k);
        state.p_last(r, c)     = p_list(k);
        active_mask(r, c)      = true;
    end

    % ----------------------------------------------------------------
    % 2. Compute per-pixel edge velocity for active pixels
    % ----------------------------------------------------------------
    % The "current time" is the timestamp of the last event in the batch
    t_now = t_list(end);

    % Initialize velocity map (seconds per pixel of edge motion)
    v_map = zeros(imgSz);

    % Get linear indices of active pixels
    [active_rows, active_cols] = find(active_mask);

    for idx = 1:numel(active_rows)
        r = active_rows(idx);
        c = active_cols(idx);
        pol = state.p_last(r, c);

        % Define observation window bounds (clamped to image edges)
        r_min = max(1, r - R);
        r_max = min(H, r + R);
        c_min = max(1, c - R);
        c_max = min(W, c + R);

        % Extract the polarity-specific timestamp patch (Eq. 4)
        if pol > 0
            patch = state.t_last_pos(r_min:r_max, c_min:c_max);
        else
            patch = state.t_last_neg(r_min:r_max, c_min:c_max);
        end

        % Flatten and sort in ascending temporal order
        timestamps = sort(patch(:), 'ascend');

        % Remove uninitialized pixels (t == 0 means never fired)
        timestamps = timestamps(timestamps > 0);

        if numel(timestamps) < seq_len
            % Not enough data to estimate velocity — fall back to using
            % whatever timestamps are available, or assign a default
            if numel(timestamps) >= 3
                % Use the full available span
                v_map(r, c) = (timestamps(end) - timestamps(1)) / ...
                    max(1, floor(numel(timestamps) / win_side));
            else
                % Too few events: assign a large default velocity
                % (conservative — long decay, preserves the pixel)
                v_map(r, c) = 0.1;  % 100 ms per pixel
            end
            continue;
        end

        % Find the position of the target pixel's timestamp in the
        % sorted sequence. The target pixel's timestamp t_i should be
        % near the middle of the sequence per the paper's formulation.
        t_i = state.t_last_any(r, c);
        [~, center_idx] = min(abs(timestamps - t_i));

        % Extract the subsequence of length seq_len centered on t_i
        % (Eq. 5): indices k to k + n*(2R+1)
        half_below = floor(seq_len / 2);
        start_idx  = center_idx - half_below;
        end_idx    = start_idx + seq_len - 1;

        % Clamp to valid range
        if start_idx < 1
            start_idx = 1;
            end_idx   = min(numel(timestamps), seq_len);
        elseif end_idx > numel(timestamps)
            end_idx   = numel(timestamps);
            start_idx = max(1, end_idx - seq_len + 1);
        end

        % Velocity (Eq. 6): time for edge to traverse n pixels
        v_map(r, c) = (timestamps(end_idx) - timestamps(start_idx)) / n;
    end

    % Ensure minimum velocity to prevent division by zero
    v_map = max(v_map, 1e-9);

    % ----------------------------------------------------------------
    % 3. Compute METS surface at t_now (Eq. 7)
    % ----------------------------------------------------------------
    % Elapsed time since last event at each pixel
    dt_map = t_now - state.t_last_any;

    % Motion-encoded decay: tau_effective = d * v(x,y)
    tau_map = d .* v_map;

    % Exponential decay
    S = exp(-dt_map ./ tau_map);

    % Hard cutoff: zero pixels where edge has moved beyond d_th pixels
    % Condition: t - t_last > d_th * v  (i.e., edge displaced > d_th px)
    cutoff_mask = dt_map > (d_th .* v_map);
    S(cutoff_mask) = 0;

    % Zero out pixels that have never fired
    S(state.t_last_any == 0) = 0;

    % Zero out pixels where velocity was not computed (not active)
    % but retain decaying values from previous events
    % (handled implicitly: v_map = 0 for inactive → tau = 0 → S = 0)

    % ----------------------------------------------------------------
    % 4. Normalize for display
    % ----------------------------------------------------------------
    % METS values are already in [0, 1] by construction (Eq. 7).
    % Unlike IEI-ATS (which is bipolar [-1,+1] and needs (S+1)/2),
    % METS is unsigned — output directly for full dynamic range.
    %
    % NOTE: The paper's figures use a warm colormap (e.g., 'hot') for
    % visualization, not raw grayscale. To replicate the paper's look:
    %   imagesc(S); colormap('hot'); axis image; caxis([0 1]);
    normalized_output_frame = S;

end