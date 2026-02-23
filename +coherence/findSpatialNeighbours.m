function [neighbor_db, D] = findSpatialNeighbours(x_valid, y_valid, t_valid, r_s, imgSz, t_interval)
% findSpatialNeighbours Find spatial-temporal neighbours within a radius.
%
% Syntax:
%   [neighbor_db, D] = findSpatialNeighbours(x_valid, y_valid, t_valid, r_s, imgSz, t_interval)
%
% Description:
%   Computes for each input point the indices and distances of neighbouring
%   points that lie within a specified spatial-temporal radius. Spatial
%   coordinates are normalized by the provided image size and temporal
%   coordinates are normalized by the provided time interval. A k-d tree
%   is used to accelerate radius queries.
%
% Inputs:
%   x_valid    - Vector of x-coordinates (pixels) for valid points (Nx1 or 1xN).
%   y_valid    - Vector of y-coordinates (pixels) for valid points (Nx1 or 1xN).
%   t_valid    - Vector of timestamps for valid points (Nx1 or 1xN).
%   r_s        - Scalar radius for neighbourhood search in normalized units.
%   imgSz      - Two-element vector [width, height] used to normalize x and y.
%                If provided as [rows, cols] the function treats imgSz(1)
%                as the x-dimension and imgSz(2) as the y-dimension.
%   t_interval - Scalar used to normalize time differences (positive).
%
% Outputs:
%   neighbor_db - Cell array (Nx1) where each cell contains indices of
%                 neighbours (including the point itself) for the
%                 corresponding input point.
%   D           - Cell array (Nx1) where each cell contains distances
%                 corresponding to the indices in neighbor_db.
%
% Notes:
%   - Inputs are cast to single precision prior to tree construction to
%     reduce memory usage.
%   - The function uses MATLAB's createns and rangesearch with a k-d tree.
%   - The function returns indices relative to the input vectors.
%
% Example:
%   % [nb, D] = findSpatialNeighbours(x, y, t, 0.05, [1024,768], 30);
%
% See also createns, rangesearch

% Validate basic inputs
assert(isvector(x_valid) && isvector(y_valid) && isvector(t_valid), ...
    'x_valid, y_valid, and t_valid must be vectors of the same length.');
assert(numel(x_valid)==numel(y_valid) && numel(x_valid)==numel(t_valid), ...
    'x_valid, y_valid, and t_valid must have the same number of elements.');
assert(isscalar(r_s) && r_s>=0, 'r_s must be a nonnegative scalar.');
assert(isvector(imgSz) && numel(imgSz)>=2, 'imgSz must be a two-element vector.');
assert(isscalar(t_interval) && t_interval>0, 't_interval must be a positive scalar.');

% Ensure column vectors
x_valid = x_valid(:);
y_valid = y_valid(:);
t_valid = t_valid(:);

% Normalize to [0,1] (use imgSz(1) for x, imgSz(2) for y)
x_norm = single(x_valid ./ imgSz(1));
y_norm = single(y_valid ./ imgSz(2));
t_norm = single((t_valid - min(t_valid)) ./ t_interval);

% Stack the data
points3D = [x_norm, y_norm, t_norm];

% Create the tree object for the rangesearch function
tree = createns(points3D, 'NSMethod', 'kdtree', 'Distance', 'euclidean', 'BucketSize', 50);

% Use MATLAB's rangesearch to find the neighbours within radius r_s
[idx, D] = rangesearch(tree, points3D, r_s);

% Ensure outputs are column cell arrays
neighbor_db = idx(:);
D = D(:);

end