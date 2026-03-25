function [last_t_map, normalized_output_frame] = ...
    timeSurface(last_t_map, x, y, t, imgSz, ts_time_constant)
% TIMESURFACE  Classical exponential-decay time surface.
%
%   [LAST_T_MAP, NORMALIZED_OUTPUT_FRAME] = TIMESURFACE(LAST_T_MAP,
%   X, Y, T, IMGSZ, TS_TIME_CONSTANT) updates the per-pixel timestamp
%   map and computes the exponential-decay time surface at t_now.
%
%   Reference:
%       Lagorce et al., "HOTS: A Hierarchy of Event-Based
%       Time-Surfaces for Pattern Recognition," IEEE T-PAMI,
%       vol. 39, no. 7, pp. 1346-1359, 2017.
%       DOI: 10.1109/TPAMI.2016.2574707
%
%   Inputs:
%     last_t_map       - [imgSz] Previous timestamp map. Initialize
%                        with -inf(imgSz) so unvisited pixels decay
%                        to zero.
%     x, y             - [N x 1] Pixel coordinates (row, col).
%     t                - [N x 1] Timestamps [s]. Scalar is broadcast.
%     imgSz            - [1 x 2] Image dimensions [nRows, nCols].
%     ts_time_constant - Scalar decay time constant tau [s].
%
%   Outputs:
%     last_t_map              - [imgSz] Updated timestamp map.
%     normalized_output_frame - [imgSz] Display surface in [0, 1].
%
%   Algorithm:
%     1. Overwrite timestamps at event locations (reset-based).
%     2. Compute surface: exp(-(t_now - t_last) / tau).
%
%   Notes:
%     - Reset-based: each event overwrites the stored timestamp.
%     - Unsigned output in [0, 1]. Direct passthrough (no rescaling).
%     - Coordinates: x = row, y = col, sub2ind(imgSz, x, y).
%
%   See also: accumulator.adaptiveGlobalDecay,
%             accumulator.localAdaptiveTimeSurface

    % ----------------------------------------------------------------
    % 0. Input validation
    % ----------------------------------------------------------------
    assert(isnumeric(ts_time_constant) ...
        && isscalar(ts_time_constant) ...
        && ts_time_constant > 0, ...
        'ts_time_constant must be a positive scalar.');

    x = x(:);
    y = y(:);
    if isscalar(t)
        t = repmat(t, numel(x), 1);
    else
        t = t(:);
        assert(numel(t) == numel(x), ...
            'Length of t must match number of coordinates.');
    end

    % ----------------------------------------------------------------
    % 1. Update timestamp map (reset-based)
    % ----------------------------------------------------------------
    linear_idx = sub2ind(imgSz, x, y);
    last_t_map(linear_idx) = t;

    % ----------------------------------------------------------------
    % 2. Compute exponential decay surface
    % ----------------------------------------------------------------
    t_now = max(t);
    decayed_surface = ...
        exp(-(t_now - last_t_map) / ts_time_constant);

    % ----------------------------------------------------------------
    % 3. Normalize for display
    % ----------------------------------------------------------------
    % Surface is unsigned [0, 1] — output directly.
    normalized_output_frame = decayed_surface;

end