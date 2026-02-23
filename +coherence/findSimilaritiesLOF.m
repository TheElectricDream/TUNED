function [similarity_score] = findSimilaritiesLOF(x_valid, y_valid, cv_valid, imgSz)
% findSimilarities Compute a simplified Local Outlier Factor-like score.
%
% similarity_score = findSimilarities(x_valid, y_valid, cv_valid, imgSz)
%
% Inputs:
%   x_valid   - vector of x coordinates (pixels)
%   y_valid   - vector of y coordinates (pixels)
%   cv_valid  - vector of a third feature (e.g., confidence or intensity)
%   imgSz     - two-element vector [height, width] of the image used to
%               normalize x and y
%
% Output:
%   similarity_score - N-by-1 vector of scores >0. Scores near 1 indicate
%                      points whose local neighborhood density is similar
%                      to their neighbors. Scores >1 indicate points that
%                      are relatively sparser (potential outliers).
%
% Notes:
% - x and y are normalized to [0,1] by dividing by imgSz. cv_valid is
%   treated as-is and cast to single precision together with x and y.
% - Uses a k-d tree for efficient nearest neighbor search and K=200
%   neighbors (including the point itself) to compute mean neighbor
%   distances. The implementation avoids divide-by-zero when duplicates
%   are present by replacing zero mean distances with eps.
% - The function returns a ratio of the mean neighborhood mean distance to
%   the point's own mean distance: values >1 suggest lower local density.
%
% Example:
%   s = findSimilarities(x, y, conf, [480, 640]);

% Normalize to [0,1]
x_norm = single(x_valid ./ imgSz(1));
y_norm = single(y_valid ./ imgSz(2));
cv_norm = single(cv_valid);

% Stack the data
points3D = [x_norm(:), y_norm(:), cv_norm(:)];

% Create the tree object for the rangesearch function
tree = createns(points3D, 'NsMethod', 'kdtree', 'Distance','euclidean', 'BucketSize', 50);

% Use MATLABs rangesearch to find the nearest neighbours
[idx, D] = knnsearch(tree, points3D, 'K', 5);

% Compute Mean Distance for every point (Inverse Density)
% D is N x K. We take the mean across the rows (neighbors).
mean_dist = mean(D, 2);

% Handle divide-by-zero if duplicates exist
mean_dist(mean_dist == 0) = eps;

% Look up the Mean Distance of the Neighbors
% idx is N x K. We want mean_dist(idx).
neighbor_mean_dists = mean_dist(idx); % Now N x K matrix

% Compute the Ratio (Simplified LOF)
% "Am I as close to my neighbors as they are to theirs?"
% If Ratio > 1, I am essentially "further out" than expected (Outlier)
% If Ratio ~ 1, I fit in well (Similar)
local_context_density = mean(neighbor_mean_dists, 2);
similarity_score = local_context_density ./ mean_dist;

end