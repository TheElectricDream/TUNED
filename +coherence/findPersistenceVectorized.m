function [closestRows, closestCols, minDists, validIndicesA] = findPersistenceVectorized(mapA, mapB, imgSz)
% findPersistenceVectorized Find nearest persistent responses between two maps
%
%   [closestRows, closestCols, minDists, validIndicesA] = ...
%       findPersistenceVectorized(mapA, mapB, imgSz)
%
%   Computes, for each active location in mapA, the closest active location
%   in mapB using a normalized 3D metric (row, column, value). The row and
%   column coordinates are normalized to the [0,1] range based on image
%   size, and the map values are used directly (assumed to be in [0,1]).
%
%   Inputs:
%     mapA  - Current map (2D array). Nonzero entries are treated as
%             targets for which nearest neighbors in mapB are sought.
%     mapB  - Reference map (2D array). Entries > 1e-6 are treated as
%             candidate persistent responses.
%     imgSz - Two-element vector [imgSz(1), imgSz(2)] specifying spatial size.
%
%   Outputs:
%     closestRows   - Row indices (linear subscripts) in mapB of the nearest
%                     neighbor for each active pixel in mapA.
%     closestCols   - Column indices in mapB corresponding to closestRows.
%     minDists      - Squared Euclidean distances (in normalized 3D space)
%                     from each active pixel in mapA to its nearest neighbor
%                     in mapB.
%     validIndicesA - Linear indices of active pixels in mapA for which a
%                     nearest neighbor was returned. Use these to assign
%                     values back into a full-sized map if needed.
%
%   Notes:
%     - If mapB contains no active entries (>1e-6) or mapA contains no
%       active entries (>0), all outputs are empty.
%     - Uses MATLAB's KD-tree (createns / knnsearch) for efficient nearest
%       neighbor queries in normalized (row, col, value) space.
%     - Returned distances are squared to allow direct comparison without
%       taking square roots when only relative distances are needed.
%
%   Example:
%     % [closestRows,closestCols,minDists,validIdx] = ...
%     %     findPersistenceVectorized(mapA, mapB, [480,640]);
%
%   Author: Alexander Crain
%   See also createns, knnsearch, ind2sub, find

% Identify the indixes in the reference map which are greater then zero so
% that we only build the map we need
idxB = find(mapB > 0);

% If there are no events in the current map, we return an empty map
if isempty(idxB)
    closestRows=[]; closestCols=[]; minDists=[]; validIndicesA=[]; return;
end

% Convert the map to linear indexing for efficiency
[b_rows, b_cols] = ind2sub([imgSz(1), imgSz(2)], idxB);
b_vals = mapB(idxB);

% Normalize Coordinates (0-1)
n_b_rows = (b_rows - 1) / (imgSz(1) - 1);
n_b_cols = (b_cols - 1) / (imgSz(2) - 1);

% The time values in both maps are assumed to be normalized already
n_b_vals = b_vals;

% Create KD-Tree 
searchSpace = [n_b_rows, n_b_cols, n_b_vals];
tree = createns(searchSpace, 'NsMethod', 'kdtree');

% Only search for points that exist in the current map
idxA = find(mapA > 0);

% If there are no points in the map, return empty
if isempty(idxA)
    closestRows=[]; closestCols=[]; minDists=[]; validIndicesA=[]; return;
end

% Convert the map to linear indexing for efficiency
[t_rows, t_cols] = ind2sub([imgSz(1), imgSz(2)], idxA);
t_vals = mapA(idxA);

% Normalize the map
n_t_rows = (t_rows - 1) / (imgSz(1) - 1);
n_t_cols = (t_cols - 1) / (imgSz(2) - 1);
n_t_vals = t_vals;

% The points from mapA correspond to the points we want to use in our
% search
queryPoints = [n_t_rows, n_t_cols, n_t_vals];

% Find THE nearest event between the current map and the reference map
[idx, d] = knnsearch(tree, queryPoints, 'K', 1);

% Calculate the output 
minDists = d.^2;  % Squared Euclidean distance
closestRows = b_rows(idx);
closestCols = b_cols(idx);
validIndicesA = idxA;  % Return linear indices to assign directly to map

end
