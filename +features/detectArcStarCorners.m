function [corner_events, sae_state] = detectArcStarCorners(...
    x, y, t, p, imgSz, sae_state, params)
%detectArcStarCorners Detect corner-events using the Arc* algorithm.
%
%   [CORNER_EVENTS, SAE_STATE] = detectArcStarCorners(X, Y, T, P, IMGSZ,
%   SAE_STATE, PARAMS) processes a batch of events in temporal order,
%   maintains the filtered Surface of Active Events (S*), and classifies
%   each event as a corner or non-corner using the Arc* algorithm.
%
%   The implementation follows Algorithm 1 from Alzugaray & Chli (2018).
%   Events are separated by polarity into independent SAEs, and the S*
%   filter (Section III-B) rejects redundant events before corner testing.
%
%   Inputs:
%       x, y    - [N x 1] vectors of pixel coordinates (1-indexed, row/col)
%                 following the IEI-ATS convention where sub2ind(imgSz,x,y)
%       t       - [N x 1] vector of timestamps [seconds]
%       p       - [N x 1] vector of polarities (0/1 or -1/+1)
%       imgSz   - [1 x 2] image dimensions [nRows, nCols] (e.g. [640, 480])
%       sae_state - Struct with persistent SAE fields:
%                   .t_last    [imgSz] - last event timestamp (any polarity)
%                   .p_last    [imgSz] - polarity of last event at pixel
%                   .t_ref_pos [imgSz] - reference timestamp, positive SAE
%                   .t_ref_neg [imgSz] - reference timestamp, negative SAE
%       params  - Struct with Arc* parameters:
%                   .kappa       - S* refractory period [s] (default 0.050)
%                   .radii       - [1 x M] circle radii (default [3, 4])
%                   .arc_bounds  - {M x 1} cell of [Lmin, Lmax] per radius
%                                  (default {[3,6], [4,8]})
%
%   Outputs:
%       corner_events - [K x 4] array of detected corners [x, y, t, p]
%       sae_state     - Updated SAE state struct
%
%   Algorithm Overview:
%       1. S* Filter (Section III-B): For each incoming event, update
%          t_last. Only update the polarity-specific reference time t_ref
%          if the event is sufficiently separated in time (t > t_last +
%          kappa) or the polarity at that pixel changed. Reject events
%          that do not update t_ref as redundant.
%
%       2. Arc* Test (Algorithm 1): For each non-redundant event, inspect
%          a circular set of elements C in the polarity-matched SAE (using
%          t_ref values). An iterative arc-growing procedure identifies
%          the contiguous arc of newest timestamps. The event is classified
%          as a corner if the arc length or its complement falls within
%          [Lmin, Lmax]. The test must pass for ALL specified radii.
%
%   Reference:
%       Alzugaray, I. and Chli, M. (2018), "Asynchronous Corner Detection
%       and Tracking for Event Cameras in Real-Time," IEEE Robotics and
%       Automation Letters, 3(4), pp. 3177-3184.
%       DOI: 10.1109/LRA.2018.2849882
%
%   See also: features.buildBresenhamCircle

    % ======================== PARAMETER SETUP ===========================
    kappa      = params.kappa;
    radii      = params.radii;
    arc_bounds = params.arc_bounds;
    n_radii    = numel(radii);

    % Precompute Bresenham circle offsets for each radius
    circle_offsets = cell(n_radii, 1);
    circle_sizes   = zeros(n_radii, 1);
    for r = 1:n_radii
        circle_offsets{r} = features.buildBresenhamCircle(radii(r));
        circle_sizes(r)   = size(circle_offsets{r}, 1);
    end

    % Ensure polarity is signed (-1 / +1)
    p = double(p(:));
    p(p == 0) = -1;

    x = double(x(:));
    y = double(y(:));
    t = double(t(:));

    N = numel(t);

    % Preallocate corner output (worst case: all events are corners)
    corner_buf = zeros(N, 4);
    n_corners  = 0;

    % Unpack SAE state for local use (avoids repeated struct access)
    t_last    = sae_state.t_last;
    p_last    = sae_state.p_last;
    t_ref_pos = sae_state.t_ref_pos;
    t_ref_neg = sae_state.t_ref_neg;

    nRows = imgSz(1);
    nCols = imgSz(2);

    % ================== PROCESS EVENTS IN TEMPORAL ORDER ================
    for i = 1:N
        ex = x(i);
        ey = y(i);
        et = t(i);
        ep = p(i);

        % ----------- S* FILTER (Section III-B) -------------------------
        % Retrieve previous state at this pixel
        prev_t = t_last(ex, ey);
        prev_p = p_last(ex, ey);

        % Check refractory condition
        polarity_changed = (prev_p ~= ep) && (prev_t > 0);
        time_separated   = (et > prev_t + kappa);

        % Always update t_last and p_last
        t_last(ex, ey) = et;
        p_last(ex, ey) = ep;

        % Only update t_ref (and proceed to detection) if filter passes.
        % The first event at a pixel (prev_t == 0) always passes.
        if ~(time_separated || polarity_changed || prev_t == 0)
            continue;  % Redundant event — skip
        end

        % Update the polarity-specific reference time and select SAE
        if ep > 0
            t_ref_pos(ex, ey) = et;
            sae_ref = t_ref_pos;
        else
            t_ref_neg(ex, ey) = et;
            sae_ref = t_ref_neg;
        end

        % ----------- ARC* CORNER TEST (Algorithm 1) --------------------
        is_corner = true;

        for r = 1:n_radii
            offsets = circle_offsets{r};
            Lmin    = arc_bounds{r}(1);
            Lmax    = arc_bounds{r}(2);
            n_circ  = circle_sizes(r);

            % --- Build circular element set C from SAE ---
            C = zeros(n_circ, 1);
            valid_circle = true;

            for ci = 1:n_circ
                cx = ex + offsets(ci, 1);
                cy = ey + offsets(ci, 2);

                % Boundary check — reject events near image border
                if cx < 1 || cx > nRows || cy < 1 || cy > nCols
                    valid_circle = false;
                    break;
                end

                C(ci) = sae_ref(cx, cy);
            end

            if ~valid_circle
                is_corner = false;
                break;
            end

            % --- Arc* iterative arc-growing (Algorithm 1) ---
            if ~arcStarTest(C, n_circ, Lmin, Lmax)
                is_corner = false;
                break;
            end

        end  % radii loop

        % --- Store corner event ---
        if is_corner
            n_corners = n_corners + 1;
            corner_buf(n_corners, :) = [ex, ey, et, ep];
        end

    end  % event loop

    % ======================== PACK OUTPUTS ==============================
    corner_events = corner_buf(1:n_corners, :);

    sae_state.t_last    = t_last;
    sae_state.p_last    = p_last;
    sae_state.t_ref_pos = t_ref_pos;
    sae_state.t_ref_neg = t_ref_neg;

end


function is_corner = arcStarTest(C, n, Lmin, Lmax)
%arcStarTest Perform the Arc* corner classification on a circular set.
%
%   This implements Algorithm 1 from Alzugaray & Chli (2018). Given a
%   circular array of SAE timestamps C, the function determines whether
%   the spatial distribution of newest timestamps forms a corner pattern.
%
%   The algorithm grows an arc of the newest contiguous elements starting
%   from the global maximum. Two supporting cursors (ECW, ECCW) traverse
%   the circle in opposite directions. The arc expands when supporting
%   elements are newer than the arc's oldest member (or during the
%   initialization phase to reach Lmin). The event is a corner if the
%   final arc length or its complement falls within [Lmin, Lmax].
%
%   Inputs:
%       C    - [n x 1] circular array of SAE timestamp values
%       n    - Number of elements in the circle
%       Lmin - Minimum arc length for a valid corner
%       Lmax - Maximum arc length for a valid corner
%
%   Output:
%       is_corner - Logical scalar

    % --- Helper: circular index wrapping ---
    next_cw  = @(idx) mod(idx, n) + 1;
    next_ccw = @(idx) mod(idx - 2, n) + 1;

    % --- Step 1: Find newest element and initialize arc ---
    [~, newest_idx] = max(C);

    % Arc boundaries (indices into C)
    arc_cw_end  = newest_idx;   % CW boundary of Anew
    arc_ccw_end = newest_idx;   % CCW boundary of Anew
    arc_len     = 1;

    % Track the oldest timestamp currently in the arc
    oldest_in_arc = C(newest_idx);

    % --- Step 2: Initialize supporting elements ---
    ecw_idx  = next_cw(newest_idx);
    eccw_idx = next_ccw(newest_idx);

    % --- Step 3: Iterate until supporting elements meet ---
    while ecw_idx ~= eccw_idx

        ecw_val  = C(ecw_idx);
        eccw_val = C(eccw_idx);

        if ecw_val >= eccw_val
            % ECW is the newer (or equal) supporting element
            if oldest_in_arc <= ecw_val || arc_len < Lmin
                % EXPAND Anew CW to include ECW (Algorithm 1, line 9)
                % ExpandUntilElement: add all elements from current CW
                % boundary to the ECW position, since ECW may have
                % advanced past the arc without prior expansion.
                [arc_cw_end, arc_len, oldest_in_arc] = ...
                    expandArcCW(C, n, arc_cw_end, ecw_idx, ...
                    arc_len, oldest_in_arc);
            end
            % Always advance ECW clockwise (Algorithm 1, line 10)
            ecw_idx = next_cw(ecw_idx);
        else
            % ECCW is the newer supporting element
            if oldest_in_arc <= eccw_val || arc_len < Lmin
                % EXPAND Anew CCW to include ECCW (Algorithm 1, line 14)
                [arc_ccw_end, arc_len, oldest_in_arc] = ...
                    expandArcCCW(C, n, arc_ccw_end, eccw_idx, ...
                    arc_len, oldest_in_arc);
            end
            % Always advance ECCW counter-clockwise (Algorithm 1, line 15)
            eccw_idx = next_ccw(eccw_idx);
        end
    end

    % --- Step 4: Classification (Algorithm 1, lines 16-19) ---
    complement_len = n - arc_len;

    is_corner = (arc_len >= Lmin && arc_len <= Lmax) || ...
                (complement_len >= Lmin && complement_len <= Lmax);
end


function [new_cw_end, new_len, new_oldest] = expandArcCW(...
    C, n, old_cw_end, target_idx, arc_len, oldest_in_arc)
%expandArcCW Expand the arc clockwise to include target_idx.
%
%   Walks CW from the current arc boundary (old_cw_end) to target_idx,
%   adding each intermediate element plus the target itself to the arc.
%   Updates arc length and the oldest timestamp in the arc.

    new_oldest = oldest_in_arc;
    idx = old_cw_end;

    while idx ~= target_idx
        idx = mod(idx, n) + 1;  % next CW
        arc_len = arc_len + 1;
        if C(idx) < new_oldest
            new_oldest = C(idx);
        end
    end

    new_cw_end = target_idx;
    new_len    = arc_len;
end


function [new_ccw_end, new_len, new_oldest] = expandArcCCW(...
    C, n, old_ccw_end, target_idx, arc_len, oldest_in_arc)
%expandArcCCW Expand the arc counter-clockwise to include target_idx.
%
%   Walks CCW from the current arc boundary (old_ccw_end) to target_idx,
%   adding each intermediate element plus the target itself to the arc.
%   Updates arc length and the oldest timestamp in the arc.

    new_oldest = oldest_in_arc;
    idx = old_ccw_end;

    while idx ~= target_idx
        idx = mod(idx - 2, n) + 1;  % next CCW
        arc_len = arc_len + 1;
        if C(idx) < new_oldest
            new_oldest = C(idx);
        end
    end

    new_ccw_end = target_idx;
    new_len     = arc_len;
end
