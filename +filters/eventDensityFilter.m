function [filter_mask, eventBuffer, n_passed, n_total] = ...
    eventDensityFilter(sorted_x, sorted_y, sorted_t, ...
    imgSz, eventBuffer, edf_params)
% EVENTDENSITYFILTER  Two-stage denoising via spatiotemporal event density.
%
%   Software implementation of the event density based denoising method
%   described in:
%
%       Y. Feng, H. Lv, H. Liu, Y. Zhang, Y. Xiao, and C. Han,
%       "Event Density Based Denoising Method for Dynamic Vision
%       Sensor," Applied Sciences, vol. 10, no. 6, p. 2024, 2020.
%       doi: 10.3390/app10062024
%
%   PRINCIPLE:
%   Real events exhibit higher spatiotemporal density than background
%   activity (BA). The filter uses a two-stage approach:
%
%     Stage 1 — Coarse filtering (random BA removal):
%       For each event, count the total number of events within an
%       L x L spatial neighbourhood over the temporal window
%       [t - dt, t]. This count is the "event density" d (Eq. 6 in
%       the paper: d = ||D||_1, the L1 norm of the density matrix).
%       If d < Psi (threshold), the event is classified as random
%       noise and rejected.
%
%     Stage 2 — Fine filtering (hot pixel removal):
%       Events that passed Stage 1 are checked for hot-pixel
%       behaviour. A 3x3 density matrix is computed from only the
%       coarse-filtered events. The centre pixel is masked out via
%       the pattern matrix P (Eq. 11), and the inner product R is
%       computed (Eq. 10). If R = 0, no other passed events exist
%       in the 8-connected neighbourhood, indicating the event is
%       an isolated hot pixel firing repeatedly at the same location.
%
%   COMPARISON WITH BAF AND STC:
%     - BAF (Delbruck 2008): Binary check — did ANY neighbour fire
%       within T? Single-event support suffices to pass.
%     - STC (Liu et al. 2015): Binary check — did the SAME cell
%       fire within dT? Single-event support suffices to pass.
%     - EDF (Feng et al. 2020): Count-based — requires MULTIPLE
%       events (>= Psi) in the L x L neighbourhood within dt.
%       This makes EDF more robust when BA density is high enough
%       that isolated noise events occasionally support each other,
%       at the cost of rejecting sparse real events at leading edges.
%
%   ADAPTATION TO FRAME-BASED PROCESSING:
%   The original algorithm processes one event at a time with a
%   lookback window. In the IEI-ATS frame-batched pipeline:
%     1. Events are re-sorted into temporal order.
%     2. A sliding-window count map is maintained: as each event
%        is processed, expired events (older than dt) are removed
%        from the count map, and the L x L neighbourhood sum is
%        read for the density check.
%     3. An event buffer carries recent events across frame
%        boundaries to provide temporal continuity at the seam.
%
%   Inputs:
%     sorted_x   - [N x 1] Row coordinates (pixel-sorted order).
%     sorted_y   - [N x 1] Column coordinates.
%     sorted_t   - [N x 1] Timestamps [s] (absolute).
%     imgSz      - [1 x 2] Sensor dimensions [nRows, nCols].
%     eventBuffer - Struct with field:
%       .x, .y, .t  - Column vectors of recent events from
%                      previous frame(s) within dt of the current
%                      frame boundary. Initialize all to empty [].
%     edf_params - Struct with fields:
%       .L       - Spatial neighbourhood side length (odd integer).
%                  Default: 5. The neighbourhood is L x L centred
%                  on the event.
%       .dt      - Temporal neighbourhood [s]. Events within
%                  [t - dt, t] are counted. Default: 5e-3 (5 ms).
%       .Psi     - Density threshold. Events with fewer than Psi
%                  neighbours in the L x L x dt volume are rejected.
%                  Default: 3.
%
%   Outputs:
%     filter_mask  - [imgSz] Binary mask. 1 = at least one event
%                    passed both stages at this pixel; 0 = rejected.
%     eventBuffer  - Updated buffer of recent events for the next
%                    frame (events within dt of the frame's end time).
%     n_passed     - Number of events passing both stages.
%     n_total      - Total events in this frame.
%
%   Notes:
%     - The density matrix D (Eq. 4) is not explicitly constructed
%       per event. Instead, a running count map is maintained and
%       the L x L subregion is summed directly. This is
%       mathematically equivalent but avoids redundant allocation.
%     - The hot pixel check (Stage 2) uses only events that passed
%       Stage 1, matching the paper's sequential pipeline.
%     - The event buffer carries at most ~dt worth of events
%       across frames. It is pruned at the start of each call.
%
%   Example:
%     % Initialize before frame loop
%     edf_params.L   = 5;
%     edf_params.dt  = 5e-3;   % 5 ms
%     edf_params.Psi = 3;
%     edf_eventBuffer.x = [];
%     edf_eventBuffer.y = [];
%     edf_eventBuffer.t = [];
%
%     % Inside frame loop
%     [filter_mask, edf_eventBuffer, n_pass, n_tot] = ...
%         filter.eventDensityFilter(sorted_x, sorted_y, ...
%         sorted_t, imgSz, edf_eventBuffer, edf_params);
%
%   See also: filter.backgroundActivityFilter,
%             filter.spatiotemporalCorrelation

    % ================================================================
    % 0. Parse parameters
    % ================================================================
    L   = edf_params.L;
    dt  = edf_params.dt;
    Psi = edf_params.Psi;

    half_L = floor(L / 2);  % Half-width of the spatial neighbourhood
    nRows  = imgSz(1);
    nCols  = imgSz(2);
    n_total = length(sorted_x);

    % ================================================================
    % 1. Re-sort current frame events into temporal order
    % ================================================================
    [temporal_t, time_order] = sort(sorted_t, 'ascend');
    temporal_x = sorted_x(time_order);
    temporal_y = sorted_y(time_order);

    % ================================================================
    % 2. Merge with event buffer from previous frame
    % ================================================================
    % The buffer contains events from the tail of the previous frame
    % that are still within dt of the current frame's start time.
    % Prepend them so they contribute to the sliding window.
    if ~isempty(eventBuffer.t)
        all_x = [eventBuffer.x; temporal_x];
        all_y = [eventBuffer.y; temporal_y];
        all_t = [eventBuffer.t; temporal_t];
        n_buffer = length(eventBuffer.t);
    else
        all_x = temporal_x;
        all_y = temporal_y;
        all_t = temporal_t;
        n_buffer = 0;
    end

    n_all = length(all_x);

    % ================================================================
    % 3. Stage 1 — Coarse filtering via sliding-window density
    % ================================================================
    % Maintain a 2D count map: how many events from the current
    % sliding window [t_i - dt, t_i] are at each pixel.
    count_map = zeros(imgSz);

    % Tail pointer: index of the oldest event still in the window
    tail = 1;

    % Pre-allocate coarse pass flags (only for current-frame events)
    coarse_pass = false(n_total, 1);

    % Seed the count map with buffer events (they are already within
    % dt of the first current-frame event, by construction)
    for i = 1:n_buffer
        count_map(all_x(i), all_y(i)) = ...
            count_map(all_x(i), all_y(i)) + 1;
    end
    % The tail starts at 1 (buffer events can expire as we advance)

    for i = (n_buffer + 1):n_all
        t_now = all_t(i);

        % --- Expire old events from the sliding window ---
        while tail <= (i - 1) && (t_now - all_t(tail)) > dt
            count_map(all_x(tail), all_y(tail)) = ...
                count_map(all_x(tail), all_y(tail)) - 1;
            tail = tail + 1;
        end

        % --- Compute event density: sum of L x L neighbourhood ---
        % Clamp the neighbourhood to array bounds
        r_lo = max(1, all_x(i) - half_L);
        r_hi = min(nRows, all_x(i) + half_L);
        c_lo = max(1, all_y(i) - half_L);
        c_hi = min(nCols, all_y(i) + half_L);

        % d = ||D||_1 (Eq. 6): sum of all events in the subregion
        d = sum(count_map(r_lo:r_hi, c_lo:c_hi), 'all');

        % --- Density threshold check (Eq. 9) ---
        frame_idx = i - n_buffer;  % Index into current-frame events
        if d >= Psi
            coarse_pass(frame_idx) = true;
        end

        % --- Add current event to the count map ---
        count_map(all_x(i), all_y(i)) = ...
            count_map(all_x(i), all_y(i)) + 1;
    end

    % ================================================================
    % 4. Stage 2 — Fine filtering (hot pixel removal)
    % ================================================================
    % Build a count map from ONLY the coarse-passed events
    coarse_count_map = zeros(imgSz);
    passed_idx = find(coarse_pass);

    for k = 1:length(passed_idx)
        j = passed_idx(k);
        coarse_count_map(temporal_x(j), temporal_y(j)) = ...
            coarse_count_map(temporal_x(j), temporal_y(j)) + 1;
    end

    % For each coarse-passed event, check 3x3 neighbourhood
    % excluding the centre pixel (Eqs. 10–11)
    % P = [1 1 1; 1 0 1; 1 1 1]
    % R = || <P, D_3x3> ||_1
    % If R = 0 → hot pixel → reject
    fine_pass = false(n_total, 1);

    for k = 1:length(passed_idx)
        j = passed_idx(k);
        ex = temporal_x(j);
        ey = temporal_y(j);

        % Extract 3x3 neighbourhood bounds
        r_lo = max(1, ex - 1);
        r_hi = min(nRows, ex + 1);
        c_lo = max(1, ey - 1);
        c_hi = min(nCols, ey + 1);

        % Sum all events in the 3x3 region
        R = sum(coarse_count_map(r_lo:r_hi, c_lo:c_hi), 'all');

        % Subtract the centre pixel's own count (the P mask zeros
        % the centre element)
        R = R - coarse_count_map(ex, ey);

        % Hot pixel check: R > 0 means at least one other passed
        % event exists in the 3x3 neighbourhood
        if R > 0
            fine_pass(j) = true;
        end
    end

    % ================================================================
    % 5. Build per-pixel binary mask from final pass results
    % ================================================================
    filter_mask = zeros(imgSz);
    n_passed = sum(fine_pass);

    final_x = temporal_x(fine_pass);
    final_y = temporal_y(fine_pass);

    if ~isempty(final_x)
        final_lin = sub2ind(imgSz, final_x, final_y);
        filter_mask(final_lin) = 1;
    end

    % ================================================================
    % 6. Update the event buffer for the next frame
    % ================================================================
    % Carry forward events from this frame that are within dt of
    % the frame's last timestamp (so they can support events at the
    % start of the next frame).
    if n_total > 0
        t_end = temporal_t(end);
        keep = temporal_t >= (t_end - dt);
        eventBuffer.x = temporal_x(keep);
        eventBuffer.y = temporal_y(keep);
        eventBuffer.t = temporal_t(keep);
    else
        eventBuffer.x = [];
        eventBuffer.y = [];
        eventBuffer.t = [];
    end

end