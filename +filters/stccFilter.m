function [filter_mask, n_passed, n_total] = ...
    stccFilter(sorted_x, sorted_y, sorted_t, sorted_p, ...
    imgSz, stcc_params)
% STCCFILTER  Space-time-content correlation filter with self-adjusting threshold.
%
%   Implementation of the STCC-Filter described in:
%
%       M. Li, Y. Huang, M. Wang, W. Li, and X. Zeng, "STCC-Filter:
%       A space-time-content correlation-based noise filter with
%       self-adjusting threshold for event camera," Signal Processing:
%       Image Communication, vol. 126, p. 117136, 2024.
%       doi: 10.1016/j.image.2024.117136
%
%   PRINCIPLE:
%   For each event, compute a Probability of Signal (POS) score by
%   convolving the N nearest events in the stream with a 2D Gaussian
%   kernel over spatial distance and temporal distance (Eq. 10).
%   Then compute a content correlation scaling factor W(p) from the
%   polarity differences between the centre event and its neighbours
%   using a 1D Gaussian kernel (Eq. 16). The discrimination threshold
%   is self-adjusted: TH' = TH / W(p). Events with POS > TH' pass;
%   the rest are rejected as noise (Eq. 18).
%
%   The content correlation addresses a key weakness of pure
%   spatiotemporal filters: under low-light conditions, BA noise
%   fires frequently enough to gain spatiotemporal support from
%   other noise events. However, noise polarities are random, while
%   real events from the same moving edge share consistent polarity.
%   By lowering the threshold for polarity-consistent events (high
%   W(p)) and raising it for polarity-inconsistent events (low W(p)),
%   the filter preserves real events that lack spatiotemporal support
%   while rejecting noise that has gained false support.
%
%   COMPARISON WITH IEI-ATS COHERENCE FILTER:
%   Both use multi-criterion scoring rather than binary pass/fail:
%     - Coherence: Three additive rules (spatial density via KD-tree,
%       temporal persistence via cross-frame KNN, IEI regularity via
%       CV). Combined score thresholded downstream.
%     - STCC: Gaussian-weighted POS score with threshold modulated
%       by polarity content correlation. Score and threshold are
%       both adaptive.
%   Key difference: Coherence uses the IEI map (inter-event interval
%   statistics) as its signal quality metric, while STCC uses polarity
%   consistency. For your EVOS dataset with signed polarity events,
%   this is a meaningful comparison axis.
%
%   ADAPTATION TO FRAME-BASED PROCESSING:
%   The paper processes N=1000 events as a batch centred on each
%   event. In the IEI-ATS frame-batched pipeline, we process all
%   events within the current frame in temporal order. For each
%   event, the N/2 preceding and N/2 following events in the
%   temporally-sorted stream form the context window. Events near
%   the start/end of the frame use asymmetric windows (all N
%   neighbours from whichever side is available). No cross-frame
%   state is needed (stateless filter).
%
%   Inputs:
%     sorted_x   - [M x 1] Row coordinates (pixel-sorted order).
%     sorted_y   - [M x 1] Column coordinates.
%     sorted_t   - [M x 1] Timestamps [s] (absolute).
%     sorted_p   - [M x 1] Signed polarity (+1 or -1).
%     imgSz      - [1 x 2] Sensor dimensions [nRows, nCols].
%     stcc_params - Struct with fields:
%       .N         - Number of neighbouring events to consider for
%                    POS and W(p) computation. Default: 1000.
%                    Larger N → more context but slower.
%       .sigma_d   - Std dev for spatial distance Gaussian [pixels].
%                    Default: 8.
%       .sigma_t   - Std dev for temporal distance Gaussian [s].
%                    Default: 5e-3 (5 ms).
%       .sigma_p   - Std dev for polarity distance Gaussian.
%                    Default: 1. Since |Δp| ∈ {0, 1}, this controls
%                    how much polarity mismatch penalises W(p).
%       .TH        - Base discrimination threshold for POS.
%                    Events pass if POS > TH / W(p).
%                    Default: 0.1 (tune to dataset).
%
%   Outputs:
%     filter_mask  - [imgSz] Binary mask. 1 = at least one event
%                    passed at this pixel; 0 = rejected.
%     n_passed     - Number of events passing the filter.
%     n_total      - Total events in this frame.
%
%   Notes:
%     - Computational cost is O(M * N) per frame, where M is the
%       number of events and N is the context window size. This is
%       the most expensive filter in the +filter/ package. For
%       frames with >10k events, consider reducing N.
%     - The polarity input must be SIGNED (+1/-1), not unsigned
%       (0/1). If your data uses 0/1 polarity, convert before
%       calling: sorted_p_signed = sorted_p * 2 - 1.
%     - The paper sets σ_p = 1 and notes that since |Δp|² ∈ {0,1},
%       the exponential for polarity mismatch reduces to a constant
%       exp(-1/(2*σ_p²)) ≈ 0.6065 when σ_p = 1. This means
%       polarity-matched neighbours contribute ~1.65x more to W(p)
%       than mismatched neighbours.
%     - The POS excludes the centre event itself (d=0, t=0) from
%       the sum to prevent self-support, consistent with the paper's
%       formulation where the sum is over the N surrounding events.
%
%   Example:
%     stcc_params.N       = 500;     % Reduce for speed
%     stcc_params.sigma_d = 8;       % pixels
%     stcc_params.sigma_t = 5e-3;    % seconds
%     stcc_params.sigma_p = 1;
%     stcc_params.TH      = 0.1;
%
%     [filter_mask, n_pass, n_tot] = ...
%         filter.stccFilter(sorted_x, sorted_y, sorted_t, ...
%         p_signed, imgSz, stcc_params);
%
%   See also: coherence.computeCoherenceMask,
%             filter.eventDensityFilter

    % ================================================================
    % 0. Parse parameters
    % ================================================================
    N       = stcc_params.N;
    sigma_d = stcc_params.sigma_d;
    sigma_t = stcc_params.sigma_t;
    sigma_p = stcc_params.sigma_p;
    TH      = stcc_params.TH;

    n_total = length(sorted_x);

    if n_total < 2
        filter_mask = zeros(imgSz);
        n_passed = 0;
        return;
    end

    % ================================================================
    % 1. Re-sort events into temporal order
    % ================================================================
    [temporal_t, time_order] = sort(sorted_t, 'ascend');
    temporal_x = sorted_x(time_order);
    temporal_y = sorted_y(time_order);
    temporal_p = sorted_p(time_order);

    % ================================================================
    % 2. Pre-compute Gaussian normalization constants
    % ================================================================
    % Spatial-temporal Gaussian (Eq. 8): 1/(2π σ_d σ_t)
    norm_st = 1 / (2 * pi * sigma_d * sigma_t);

    % Two pre-computed denominators for the exponent
    inv_2sigma_d2 = 1 / (2 * sigma_d^2);
    inv_2sigma_t2 = 1 / (2 * sigma_t^2);

    % Content Gaussian (Eq. 14): 1/(√(2π) σ_p)
    norm_p = 1 / (sqrt(2 * pi) * sigma_p);
    inv_2sigma_p2 = 1 / (2 * sigma_p^2);

    % Pre-compute the two possible G(p) values since |Δp|² ∈ {0, 1}
    % When Δp = 0 (same polarity): G(p) = norm_p * exp(0) = norm_p
    % When |Δp| = 1 (diff polarity): G(p) = norm_p * exp(-1/(2σ_p²))
    % For unsigned polarity {0,1}: |Δp| ∈ {0,1}
    % For signed polarity {-1,+1}: |Δp| ∈ {0,2}, so (Δp)² ∈ {0,4}
    % The paper uses p ∈ {0,1}, but we support signed polarity.
    % We compute Δp_i = |p_i - p_0| for each pair and use it directly.

    % ================================================================
    % 3. Per-event POS and W(p) computation
    % ================================================================
    half_N = floor(N / 2);
    pass_flag = false(n_total, 1);
    n_passed = 0;

    for i = 1:n_total
        % Centre event
        x0 = temporal_x(i);
        y0 = temporal_y(i);
        t0 = temporal_t(i);
        p0 = temporal_p(i);

        % Determine the context window: N/2 before, N/2 after
        j_start = max(1, i - half_N);
        j_end   = min(n_total, i + half_N);

        % If one side is truncated, extend the other
        if (i - j_start) < half_N
            j_end = min(n_total, j_start + N);
        elseif (j_end - i) < half_N
            j_start = max(1, j_end - N);
        end

        % --- Vectorised POS computation (Eq. 10) ---
        % Extract neighbour coordinates (excluding self)
        idx = [j_start:(i-1), (i+1):j_end]';

        if isempty(idx)
            % Isolated event — no neighbours
            continue;
        end

        dx = double(temporal_x(idx)) - double(x0);
        dy = double(temporal_y(idx)) - double(y0);
        dt = double(temporal_t(idx)) - double(t0);

        % Spatial distance (Eq. 3): Euclidean
        dd_sq = dx.^2 + dy.^2;

        % Temporal distance squared (Eq. 4)
        dt_sq = dt.^2;

        % Gaussian-weighted POS (Eq. 10)
        gauss_st = norm_st .* exp(-(dd_sq .* inv_2sigma_d2 + ...
                                     dt_sq .* inv_2sigma_t2));
        POS = sum(gauss_st);

        % --- Vectorised W(p) computation (Eq. 16) ---
        dp = abs(double(temporal_p(idx)) - double(p0));
        dp_sq = dp.^2;

        gauss_p = norm_p .* exp(-dp_sq .* inv_2sigma_p2);
        Wp = sum(gauss_p);

        % --- Self-adjusting threshold (Eq. 17–18) ---
        % TH' = TH / W(p)
        % Event passes if POS > TH'
        % Guard against W(p) = 0 (all neighbours have extreme
        % polarity mismatch — effectively infinite threshold)
        if Wp > 0
            TH_adj = TH / Wp;
        else
            TH_adj = Inf;  % Cannot pass
        end

        if POS > TH_adj
            pass_flag(i) = true;
            n_passed = n_passed + 1;
        end
    end

    % ================================================================
    % 4. Build per-pixel binary mask
    % ================================================================
    filter_mask = zeros(imgSz);

    passed_x = temporal_x(pass_flag);
    passed_y = temporal_y(pass_flag);

    if ~isempty(passed_x)
        passed_lin = sub2ind(imgSz, passed_x, passed_y);
        filter_mask(passed_lin) = 1;
    end

end