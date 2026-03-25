function [filter_mask, lastTimesMap, n_passed, n_total] = ...
    backgroundActivityFilter(sorted_x, sorted_y, sorted_t, ...
    imgSz, lastTimesMap, baf_params)
% BACKGROUNDACTIVITYFILTER  Neighbour-support BA noise filter.
%
%   Software implementation of the Background Activity Filter
%   described in:
%
%       T. Delbruck, "Frame-free dynamic digital vision," Proc. Intl.
%       Symp. on Secure-Life Electronics, Advanced Electronics for
%       Quality Life and Society, University of Tokyo, 2008, pp. 21–26.
%
%   Also implemented in the jAER open-source project as
%   BackgroundActivityFilter (see jAER [1]).
%
%   PRINCIPLE:
%   Background activity (BA) events are temporally uncorrelated with
%   their spatial neighbourhood, while real events from moving edges
%   are spatiotemporally correlated. The filter exploits this by
%   maintaining a per-pixel timestamp map and using a neighbour-write
%   strategy:
%
%     1. CHECK: Read the stored timestamp at the current event's
%        pixel (x, y). If the time difference is within the support
%        window T, a neighbour must have recently written there →
%        PASS the event. Otherwise → REJECT.
%
%     2. WRITE: Store the current event's timestamp into all 8
%        neighbouring pixel locations (not the event's own pixel).
%        This is the key trick from Delbruck (2008, Section IV):
%        "This implementation avoids iteration and branching over
%        all neighboring pixels by simply storing an event's
%        timestamp in all neighbors. Then only a single conditional
%        branch is necessary."
%
%   This means the timestamp stored at any pixel was actually
%   written by one of its 8 neighbours. The single read-and-compare
%   at the current pixel implicitly checks whether *any* neighbour
%   had recent activity.
%
%   COMPARISON WITH STC FILTER (Liu et al. 2015):
%   The STC filter checks the inter-spike interval at the *same*
%   subsampled cell. Spatial support comes from address subsampling
%   (grouping NxN pixel blocks). The BA filter instead provides
%   spatial support through the 8-connected neighbour write pattern,
%   operating at full pixel resolution. The BA filter is simpler
%   (no subsampling logic) but has a fixed 3x3 spatial footprint,
%   while the STC filter's spatial neighbourhood is programmable.
%
%   ADAPTATION TO FRAME-BASED PROCESSING:
%   As with the STC filter implementation, events are re-sorted
%   into temporal order within each frame window to preserve causal
%   ISI semantics. The lastTimesMap persists across frames.
%
%   Inputs:
%     sorted_x   - [N x 1] Row coordinates of events in current frame,
%                  sorted by linear pixel index (IEI-ATS convention).
%     sorted_y   - [N x 1] Column coordinates of events.
%     sorted_t   - [N x 1] Timestamps [s] (absolute).
%     imgSz      - [1 x 2] Sensor dimensions [nRows, nCols].
%     lastTimesMap - [imgSz] Persistent timestamp map at full pixel
%                  resolution. Initialize to -Inf on first call.
%     baf_params - Struct with fields:
%       .T       - Support time window [s]. An event passes if the
%                  timestamp stored at its pixel (written by a
%                  neighbour) is within T of the current time.
%                  Typical values: 1e-3 to 100e-3 (1–100 ms).
%
%   Outputs:
%     filter_mask   - [imgSz] Binary mask at full pixel resolution.
%                     1 = at least one correlated event passed at
%                     this pixel; 0 = no events passed.
%     lastTimesMap  - [imgSz] Updated timestamp map for next frame.
%     n_passed      - Number of events that passed the filter.
%     n_total       - Total number of events in this frame.
%
%   Notes:
%     - Boundary handling: neighbour writes that would fall outside
%       the sensor array are silently skipped (no wrap-around).
%     - The first event at any location will always fail unless a
%       neighbour has already fired within T. This is correct
%       behaviour: an isolated first event has no spatial support.
%     - Unlike the STC filter, the current event does NOT write to
%       its own pixel. It only writes to its 8 neighbours. This
%       ensures the check at (x,y) always reflects neighbour
%       activity, not self-activity.
%     - Performance: the per-event loop writes 8 neighbour entries
%       and does 1 read + 1 compare. This is O(N) per frame with
%       a small constant factor, matching Delbruck's reported
%       0.1 µs/event on 2005 hardware.
%
%   References:
%     [1] T. Delbruck, "jAER open source project," 2007.
%         Available: http://jaer.wiki.sourceforge.net
%
%   Example:
%     % Initialize before frame loop
%     baf_params.T = 10e-3;  % 10 ms support window
%     baf_lastTimesMap = -Inf(imgSz);
%
%     % Inside frame loop
%     [filter_mask, baf_lastTimesMap, n_pass, n_tot] = ...
%         filter.backgroundActivityFilter(sorted_x, sorted_y, ...
%         sorted_t, imgSz, baf_lastTimesMap, baf_params);
%
%   See also: filter.spatiotemporalCorrelation,
%             coherence.computeCoherenceMask

    % ================================================================
    % 0. Parse parameters
    % ================================================================
    T_support = baf_params.T;

    nRows = imgSz(1);
    nCols = imgSz(2);
    n_total = length(sorted_x);

    % ================================================================
    % 1. Re-sort events into temporal order
    % ================================================================
    % The input events are pixel-sorted for the IEI-ATS statistics
    % pipeline. The BA filter requires causal processing: each
    % event's neighbour-write must precede subsequent events'
    % read-check at those neighbour locations.
    [temporal_t, time_order] = sort(sorted_t, 'ascend');
    temporal_x = sorted_x(time_order);
    temporal_y = sorted_y(time_order);

    % ================================================================
    % 2. Pre-compute the 8-connected neighbour offsets
    % ================================================================
    % Relative (row, col) offsets for the 8-connected neighbourhood.
    % The event writes its timestamp to all 8 neighbours but NOT to
    % its own pixel — this is essential to the Delbruck trick.
    dx = [-1, -1, -1,  0, 0,  1, 1, 1];
    dy = [-1,  0,  1, -1, 1, -1, 0, 1];

    % ================================================================
    % 3. Sequential event processing
    % ================================================================
    pass_flag = false(n_total, 1);
    n_passed  = 0;

    for i = 1:n_total
        ex = temporal_x(i);
        ey = temporal_y(i);
        t_now = temporal_t(i);

        % ----------------------------------------------------------
        % Step 2 (Delbruck Section IV): CHECK
        % Read the timestamp at this event's own pixel. If it is
        % within T of the current time, a neighbour wrote there
        % recently → this event is spatiotemporally correlated.
        % ----------------------------------------------------------
        dt = t_now - lastTimesMap(ex, ey);

        if dt <= T_support
            pass_flag(i) = true;
            n_passed = n_passed + 1;
        end

        % ----------------------------------------------------------
        % Step 1 (Delbruck Section IV): WRITE
        % Store the event's timestamp in all 8 neighbouring pixel
        % locations. Boundary-safe: skip writes outside the array.
        % ----------------------------------------------------------
        for k = 1:8
            nx = ex + dx(k);
            ny = ey + dy(k);

            if nx >= 1 && nx <= nRows && ny >= 1 && ny <= nCols
                lastTimesMap(nx, ny) = t_now;
            end
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