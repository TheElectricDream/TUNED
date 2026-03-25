function [closestRows, closestCols, minDists, validIndicesA] = ...
    findPersistenceVectorized(mapA, mapB, imgSz)
% FINDPERSISTENCEVECTORIZED  Cross-frame persistence via KNN search.
%
%   [CLOSESTROWS, CLOSESTCOLS, MINDISTS, VALIDINDICESA] =
%   FINDPERSISTENCEVECTORIZED(MAPA, MAPB, IMGSZ) finds, for each
%   active pixel in mapA, the nearest active pixel in mapB using a
%   normalized 3D metric over (row, col, value) space. This measures
%   whether active regions persist across consecutive frames.
%
%   Inputs:
%     mapA  - [imgSz] Current frame map. Nonzero entries are query
%             points.
%     mapB  - [imgSz] Previous frame map. Entries > 1e-6 are
%             candidate persistent responses.
%     imgSz - [1 x 2] Image dimensions [nRows, nCols].
%
%   Outputs:
%     closestRows   - [M x 1] Row indices in mapB of nearest
%                     neighbour for each query point.
%     closestCols   - [M x 1] Column indices in mapB.
%     minDists      - [M x 1] Squared Euclidean distances in
%                     normalized (row, col, value) space.
%     validIndicesA - [M x 1] Linear indices of active pixels in
%                     mapA for which a neighbour was found.
%
%   Algorithm:
%     1. Extract active pixels from mapB (reference) and normalize
%        coordinates to [0, 1].
%     2. Build KD-tree over the reference points.
%     3. Extract active pixels from mapA (query) and normalize.
%     4. K=1 nearest neighbour search from query to reference.
%     5. Return squared distances for direct threshold comparison.
%
%   Notes:
%     - Returns empty outputs if either map has no active pixels.
%     - Squared distances avoid unnecessary sqrt operations when
%       only relative ordering or threshold comparison is needed.
%     - Coordinates: row/col normalized by imgSz; values assumed
%       to already be in [0, 1] from upstream log-normalization.
%
%   See also: coherence.computeCoherenceMask, createns, knnsearch

    % ----------------------------------------------------------------
    % 0. Extract reference points from mapB
    % ----------------------------------------------------------------
    idxB = find(mapB > 1e-6);

    if isempty(idxB)
        closestRows = [];
        closestCols = [];
        minDists = [];
        validIndicesA = [];
        return;
    end

    [b_rows, b_cols] = ind2sub(imgSz, idxB);
    b_vals = mapB(idxB);

    % Normalize to [0, 1]
    n_b_rows = (b_rows - 1) / (imgSz(1) - 1);
    n_b_cols = (b_cols - 1) / (imgSz(2) - 1);

    % ----------------------------------------------------------------
    % 1. Build KD-tree over reference points
    % ----------------------------------------------------------------
    searchSpace = [n_b_rows, n_b_cols, b_vals];
    tree = createns(searchSpace, 'NsMethod', 'kdtree');

    % ----------------------------------------------------------------
    % 2. Extract query points from mapA
    % ----------------------------------------------------------------
    idxA = find(mapA > 0);

    if isempty(idxA)
        closestRows = [];
        closestCols = [];
        minDists = [];
        validIndicesA = [];
        return;
    end

    [t_rows, t_cols] = ind2sub(imgSz, idxA);
    t_vals = mapA(idxA);

    n_t_rows = (t_rows - 1) / (imgSz(1) - 1);
    n_t_cols = (t_cols - 1) / (imgSz(2) - 1);

    queryPoints = [n_t_rows, n_t_cols, t_vals];

    % ----------------------------------------------------------------
    % 3. K=1 nearest neighbour search
    % ----------------------------------------------------------------
    [idx, d] = knnsearch(tree, queryPoints, 'K', 1);

    minDists = d .^ 2;
    closestRows = b_rows(idx);
    closestCols = b_cols(idx);
    validIndicesA = idxA;

end