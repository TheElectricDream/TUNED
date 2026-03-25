function [filter_mask, lastTimesMap, n_passed, n_total] = ...
    spatiotemporalCorrelation(sorted_x, sorted_y, sorted_t, ...
    imgSz, lastTimesMap, stc_params)
% SPATIOTEMPORALCORRELATION  Background activity filter via ISI thresholding.
%
%   Software implementation of the spatiotemporal correlation filter
%   described in:
%
%       H. Liu, C. Brandli, C. Li, S.-C. Liu, and T. Delbruck,
%       "Design of a Spatiotemporal Correlation Filter for Event-based
%       Sensors," Proc. IEEE Int. Symp. Circuits Syst. (ISCAS), 2015,
%       pp. 722–725. doi: 10.1109/ISCAS.2015.7168735
%
%   The filter exploits the key difference between real activity and
%   background activity (BA): real events exhibit spatiotemporal
%   correlation with their neighbourhood, while BA events are
%   temporally uncorrelated. For each incoming event, the inter-spike
%   interval (ISI) — time since the last event at the same subsampled
%   cell — is compared against a programmable time window dT. Events
%   with ISI < dT are passed as correlated; events with ISI >= dT
%   are rejected as likely BA.
%
%   Spatial subsampling groups blocks of sensor pixels onto the same
%   filter cell, providing the spatial correlation neighbourhood.
%   This mirrors the hardware's address-LSB-masking scheme (Fig. 2
%   in the paper).
%
%   ADAPTATION TO FRAME-BASED PROCESSING:
%   The original filter is event-driven (processes one event at a
%   time). This implementation adapts it to the IEI-ATS frame-based
%   pipeline by:
%     1. Re-sorting the pixel-sorted events back into temporal order
%        within the current frame window.
%     2. Sequentially processing each event against the persistent
%        lastTimesMap (which carries state across frames).
%     3. Building a per-pixel binary mask from the pass/fail results.
%
%   Inputs:
%     sorted_x   - [N x 1] Row coordinates of events in current frame,
%                  sorted by linear pixel index (as used by the IEI-ATS
%                  statistics pipeline). Will be re-sorted by time.
%     sorted_y   - [N x 1] Column coordinates of events.
%     sorted_t   - [N x 1] Timestamps [s] (absolute, from recording
%                  start), sorted to match sorted_x/y.
%     imgSz      - [1 x 2] Sensor dimensions [nRows, nCols].
%     lastTimesMap - [cellRows x cellCols] Persistent timestamp map at
%                  cell resolution. Initialize to -Inf on first call.
%                  Carries across frames to maintain ISI continuity.
%     stc_params - Struct with fields:
%       .dT              - Correlation time window [s]. Events with
%                          ISI < dT pass the filter. Equivalent to
%                          C*(Vrs-Vth)/I1 in the hardware. Typical
%                          values: 1e-3 to 10e-3 (1–10 ms).
%       .subsample_rate  - Address LSB masking rate. Controls the
%                          spatial neighbourhood size:
%                            0 → 1x1  (per-pixel, no subsampling)
%                            1 → 2x2  blocks
%                            2 → 4x4  blocks
%                            3 → 8x8  blocks
%                          Default: 1 (2x2, matching paper's example).
%
%   Outputs:
%     filter_mask   - [imgSz] Binary mask at full pixel resolution.
%                     1 = at least one correlated event passed at
%                     this pixel; 0 = no events passed (or no events).
%     lastTimesMap  - [cellRows x cellCols] Updated timestamp map
%                     for the next frame.
%     n_passed      - Scalar. Number of events that passed the filter.
%     n_total       - Scalar. Total number of events in this frame.
%
%   Notes:
%     - The first event at any cell will always be rejected (no
%       prior timestamp to compare against). This matches the
%       hardware behaviour shown in Fig. 5, where e1 is always
%       filtered away.
%     - All events (pass and fail) update the lastTimesMap. This
%       matches the hardware reset cycle: Vcap is reset to Vrs on
%       every Ack, regardless of the Pass signal.
%     - The filter_mask is binary (not continuous like the coherence
%       score). Downstream smoothing in main.m can be applied if
%       desired.
%     - For evaluation against the coherence filter, the filter_mask
%       is the direct comparison target: both produce a 2D mask that
%       gates the accumulator.
%
%   Example:
%     % Initialize before frame loop
%     stc_params.dT = 5e-3;            % 5 ms correlation window
%     stc_params.subsample_rate = 1;   % 2x2 spatial subsampling
%     block_size = 2^stc_params.subsample_rate;
%     cellSz = ceil(imgSz ./ block_size);
%     lastTimesMap = -Inf(cellSz);
%
%     % Inside frame loop (after pixel-sorting)
%     [filter_mask, lastTimesMap, n_pass, n_tot] = ...
%         filter.spatiotemporalCorrelation(sorted_x, sorted_y, ...
%         sorted_t, imgSz, lastTimesMap, stc_params);
%
%   See also: coherence.computeCoherenceMask, main

    % ================================================================
    % 0. Parse parameters
    % ================================================================
    dT              = stc_params.dT;
    subsample_rate  = stc_params.subsample_rate;
    block_size      = 2^subsample_rate;

    n_total = length(sorted_x);

    % ================================================================
    % 1. Re-sort events into temporal order
    % ================================================================
    % The input events are sorted by linear pixel index (for the
    % statistics and coherence pipeline). The STC filter requires
    % causal temporal processing — each event's ISI depends on all
    % preceding events at its cell, including those earlier in this
    % same frame window.
    [temporal_t, time_order] = sort(sorted_t, 'ascend');
    temporal_x = sorted_x(time_order);
    temporal_y = sorted_y(time_order);

    % ================================================================
    % 2. Compute subsampled cell indices for all events (vectorized)
    % ================================================================
    % Mirrors the hardware's address-LSB-masking: floor-divide
    % coordinates by block_size, then add 1 for MATLAB 1-indexing.
    cell_x = floor((temporal_x - 1) / block_size) + 1;
    cell_y = floor((temporal_y - 1) / block_size) + 1;

    % ================================================================
    % 3. Sequential ISI check (event-by-event)
    % ================================================================
    % Pre-allocate per-event pass flag
    pass_flag = false(n_total, 1);
    n_passed  = 0;

    for i = 1:n_total
        cx = cell_x(i);
        cy = cell_y(i);
        t_now = temporal_t(i);

        % Compute ISI: time since last event at this cell
        isi = t_now - lastTimesMap(cx, cy);

        % Correlation check (Fig. 5 in paper: e2 passes, e1/e3 fail)
        if isi < dT
            pass_flag(i) = true;
            n_passed = n_passed + 1;
        end

        % Reset: update cell timestamp regardless of pass/fail
        % (matches hardware Ack cycle — Vcap resets on every event)
        lastTimesMap(cx, cy) = t_now;
    end

    % ================================================================
    % 4. Build per-pixel binary mask from pass results
    % ================================================================
    % A pixel is marked as 1 if *any* event at that pixel passed
    % the correlation check within this frame.
    filter_mask = zeros(imgSz);

    % Extract coordinates of passed events
    passed_x = temporal_x(pass_flag);
    passed_y = temporal_y(pass_flag);

    if ~isempty(passed_x)
        % Use linear indexing for vectorized assignment
        passed_lin = sub2ind(imgSz, passed_x, passed_y);

        % Mark all passed pixel locations
        filter_mask(passed_lin) = 1;
    end

end