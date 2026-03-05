function [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, current_idx)
% SLICETOVALIDRANGE  Extract and filter events within a time window.
%
%   [CURRENT_IDX, X_VALID, Y_VALID, T_VALID, P_VALID] =
%   SLICETOVALIDRANGE(T_RANGE_N, XK, YK, TK, PK, IMGSZ,
%   CURRENT_IDX) extracts events with timestamps up to t_range_n
%   starting from current_idx, then filters to valid spatial bounds.
%
%   Inputs:
%     t_range_n   - Scalar upper time bound for this window [s].
%     xk, yk      - [M x 1] Full event coordinate vectors (row, col).
%     tk          - [M x 1] Full timestamp vector [s], sorted ascending.
%     pk          - [M x 1] Full polarity vector.
%     imgSz       - [1 x 2] Image dimensions [nRows, nCols].
%     current_idx - Scalar start index into the full vectors.
%
%   Outputs:
%     current_idx - Updated index (points to first event beyond
%                   this window, for the next call).
%     x_valid, y_valid, t_valid, p_valid - Filtered event vectors
%                   within both the time window and spatial bounds.
%
%   Algorithm:
%     1. Linear scan from current_idx to find events <= t_range_n.
%     2. Slice the event vectors.
%     3. Filter to valid spatial bounds [1, imgSz].
%
%   Notes:
%     - Assumes tk is sorted ascending — uses linear scan, not
%       binary search, for O(N_window) efficiency per call.
%     - Returns empty arrays if no events fall in the window.
%     - Coordinates: x = row, y = col.
%
%   See also: main

    % ----------------------------------------------------------------
    % 1. Linear scan to find window bounds
    % ----------------------------------------------------------------
    start_node = current_idx;
    while current_idx <= length(tk) && tk(current_idx) <= t_range_n
        current_idx = current_idx + 1;
    end
    end_node = current_idx - 1;

    % ----------------------------------------------------------------
    % 2. Slice and filter
    % ----------------------------------------------------------------
    if end_node >= start_node
        x_slice = xk(start_node:end_node);
        y_slice = yk(start_node:end_node);
        t_slice = tk(start_node:end_node);
        p_slice = pk(start_node:end_node);

        % Spatial bounds filter
        valid_mask = x_slice >= 1 & x_slice <= imgSz(1) & ...
                     y_slice >= 1 & y_slice <= imgSz(2);

        x_valid = x_slice(valid_mask);
        y_valid = y_slice(valid_mask);
        t_valid = t_slice(valid_mask);
        p_valid = p_slice(valid_mask);
    else
        x_valid = [];
        y_valid = [];
        t_valid = [];
        p_valid = [];
    end

end