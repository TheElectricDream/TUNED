function [neighbor_db, D] = findSpatialNeighbours(...
    x_valid, y_valid, t_valid, r_s, imgSz, t_interval)
% FINDSPATIALNEIGHBOURS  KD-tree radius search for spatial density.
%
%   [NEIGHBOR_DB, D] = FINDSPATIALNEIGHBOURS(X_VALID, Y_VALID,
%   T_VALID, R_S, IMGSZ, T_INTERVAL) computes for each input event
%   the indices and distances of neighbouring events within a
%   specified radius in normalized spatio-temporal (x, y, t) space.
%
%   Inputs:
%     x_valid    - [N x 1] Row coordinates (pixels).
%     y_valid    - [N x 1] Column coordinates (pixels).
%     t_valid    - [N x 1] Timestamps [s].
%     r_s        - Scalar search radius in normalized units.
%     imgSz      - [1 x 2] Image dimensions [nRows, nCols].
%                  Used to normalize spatial coordinates to [0, 1].
%     t_interval - Scalar frame interval [s]. Used to normalize
%                  temporal coordinates.
%
%   Outputs:
%     neighbor_db - {N x 1} Cell array of neighbour indices per event
%                   (including the event itself).
%     D           - {N x 1} Cell array of Euclidean distances
%                   corresponding to neighbor_db.
%
%   Algorithm:
%     1. Normalize coordinates: x/imgSz(1), y/imgSz(2),
%        (t - t_min)/t_interval.
%     2. Build KD-tree over the 3D normalized points.
%     3. Execute rangesearch with radius r_s.
%
%   Notes:
%     - Inputs are cast to single precision to reduce memory.
%     - Coordinates: x = row, y = col.
%
%   See also: coherence.computeCoherenceMask, createns, rangesearch

    % ----------------------------------------------------------------
    % 0. Input validation
    % ----------------------------------------------------------------
    assert(isvector(x_valid) && isvector(y_valid) ...
        && isvector(t_valid), ...
        'x_valid, y_valid, t_valid must be vectors.');
    assert(numel(x_valid) == numel(y_valid) ...
        && numel(x_valid) == numel(t_valid), ...
        'x_valid, y_valid, t_valid must have the same length.');
    assert(isscalar(r_s) && r_s >= 0, ...
        'r_s must be a nonnegative scalar.');
    assert(isvector(imgSz) && numel(imgSz) >= 2, ...
        'imgSz must be a two-element vector.');
    assert(isscalar(t_interval) && t_interval > 0, ...
        't_interval must be a positive scalar.');

    % ----------------------------------------------------------------
    % 1. Normalize to [0, 1] and build KD-tree
    % ----------------------------------------------------------------
    x_valid = x_valid(:);
    y_valid = y_valid(:);
    t_valid = t_valid(:);

    x_norm = single(x_valid ./ imgSz(1));
    y_norm = single(y_valid ./ imgSz(2));
    t_norm = single((t_valid - min(t_valid)) ./ t_interval);

    points3D = [x_norm, y_norm, t_norm];

    tree = createns(points3D, 'NSMethod', 'kdtree', ...
        'Distance', 'euclidean', 'BucketSize', 50);

    % ----------------------------------------------------------------
    % 2. Radius search
    % ----------------------------------------------------------------
    [idx, D] = rangesearch(tree, points3D, r_s);

    neighbor_db = idx(:);
    D = D(:);

end