function [last_t_map, normalized_output_frame] = timeSurface(last_t_map, x, y, t, imgSz, ts_time_constant)
%timeSurface Update timestamp map and compute exponential decayed surface.
%   [LAST_T_MAP, NORMALIZED_OUTPUT_FRAME] = TIMESURFACE(LAST_T_MAP, X, Y, T, IMGSZ, TS_TIME_CONSTANT)
%   updates the timestamp map LAST_T_MAP at coordinates (X,Y) with timestamps T
%   and returns the normalized exponential decayed surface NORMALIZED_OUTPUT_FRAME computed using
%   time constant TS_TIME_CONSTANT. IMGSZ is the image size used for indexing.
%
%   Inputs:
%     LAST_T_MAP        - matrix of previous timestamps sized IMGSZ
%     X, Y              - vectors of row and column coordinates (same length)
%     T                 - vector of timestamps corresponding to (X,Y)
%     IMGSZ             - two-element vector [nrows, ncols]
%     TS_TIME_CONSTANT  - positive scalar time constant for exponential decay
%
%   Outputs:
%     LAST_T_MAP        - updated timestamp map
%     DECAYED_SURFACE   - matrix of same size as LAST_T_MAP containing
%                         exp(-(t_now - LAST_T_MAP)/TS_TIME_CONSTANT)
%
%   Notes:
%     - Coordinates X and Y are expected to be within the image bounds.
%     - The function uses linear indexing via SUB2IND for assignment.

    % Validate inputs minimally
    assert(isnumeric(ts_time_constant) && isscalar(ts_time_constant) && ts_time_constant > 0, ...
        'ts_time_constant must be a positive scalar.');

    % Ensure column vectors for indexing
    x = x(:); y = y(:);
    if isscalar(t)
        t = repmat(t, numel(x), 1);
    else
        t = t(:);
        assert(numel(t) == numel(x), 'Length of t must match number of coordinates.');
    end

    % Update the Global State Maps
    linear_idx = sub2ind(imgSz, x, y);

    % Update timestamps (latest event overwrites previous ones)
    last_t_map(linear_idx) = t;

    % Define 'now'. Usually, this is the timestamp of the very last event processed.
    t_now = max(t);

    % Calculate the actual surface
    decayed_surface = exp(-(t_now - last_t_map) / ts_time_constant);
    
    % Normalize the surface
    normalized_output_frame = (decayed_surface + 1) / 2;
end