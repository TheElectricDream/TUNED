function [matches, scores] = matchCVDescriptors(desc_curr, corners_curr, ...
    desc_prev, corners_prev, varargin)
%matchCVDescriptors Match CV-patch descriptors between two frames.
%
%   [MATCHES, SCORES] = matchCVDescriptors(DESC_CURR, CORNERS_CURR, ...
%       DESC_PREV, CORNERS_PREV)
%   [MATCHES, SCORES] = matchCVDescriptors(..., Name, Value)
%
%   Finds correspondences between Harris-CV corners detected in two
%   consecutive frames by matching their local CV-patch descriptors.
%   Matching uses Normalised Cross-Correlation (NCC) — since descriptors
%   are zero-mean, unit-norm, NCC = dot product — with three rejection
%   criteria:
%
%     1. NCC threshold: reject matches below a minimum similarity.
%     2. Lowe's ratio test: reject matches where the best and second-best
%        NCC scores are too similar (ambiguous matches).
%     3. Spatial gate: reject matches where the corner moved more than
%        a maximum distance between frames (physically implausible
%        for the given frame rate).
%
%   Inputs:
%     desc_curr    - [M x D] current-frame descriptors (from
%                    extractCVDescriptors).
%     corners_curr - [M x 2] current-frame corner positions [row, col].
%     desc_prev    - [N x D] previous-frame descriptors.
%     corners_prev - [N x 2] previous-frame corner positions [row, col].
%
%   Name-Value Parameters:
%     'ncc_threshold'  - Minimum NCC score to accept a match.
%                        Range [-1, 1]. Default: 0.7
%     'ratio_test'     - Lowe's ratio threshold. A match is accepted
%                        only if best_score / second_best < ratio_test.
%                        Set to 1.0 to disable. Default: 0.8
%                        Reference: Lowe (2004), IJCV 60(2), pp. 91–110.
%     'max_distance'   - Maximum spatial displacement [pixels] between
%                        matched corners. Rejects physically implausible
%                        correspondences. Default: 50
%     'mutual'         - If true, enforce mutual (symmetric) matching:
%                        a match (i,j) is kept only if j's best match
%                        in the current frame is also i.
%                        Default: true
%
%   Outputs:
%     matches - [P x 2] array of matched index pairs:
%               matches(k,1) = index into desc_curr / corners_curr
%               matches(k,2) = index into desc_prev / corners_prev
%               Sorted by descending NCC score.
%     scores  - [P x 1] NCC scores for each accepted match.
%
%   Example:
%     % Frame N:
%     [desc_N, vc_N] = features.extractCVDescriptors(cv_map_N, corners_N);
%     % Frame N-1 (stored from previous iteration):
%     [matches, scores] = features.matchCVDescriptors( ...
%         desc_N, vc_N, desc_prev, vc_prev, ...
%         'ncc_threshold', 0.7, 'ratio_test', 0.8, 'max_distance', 40);
%     % matches(:,1) indexes into vc_N, matches(:,2) into vc_prev
%
%   Notes:
%     - The spatial gate is essential for event camera data where the
%       CV field can have repeated patterns (e.g., two similar corners
%       on opposite ends of a spacecraft). Without it, a corner on the
%       leading edge could falsely match to the trailing edge.
%     - Lowe's ratio test (IJCV 2004) is the standard ambiguity
%       rejection mechanism for local feature matching. A threshold of
%       0.8 rejects ~90% of false matches while retaining ~95% of
%       correct matches (Lowe's original finding on SIFT).
%
%   References:
%     Lowe, D.G. (2004), "Distinctive Image Features from Scale-Invariant
%       Keypoints," IJCV, 60(2), pp. 91–110.
%
%   See also: features.extractCVDescriptors, features.detectHarrisCV

    % ----------------------------------------------------------------
    % 0. Parse inputs
    % ----------------------------------------------------------------
    ip = inputParser;
    addRequired(ip, 'desc_curr',    @(x) isnumeric(x) && ismatrix(x));
    addRequired(ip, 'corners_curr', @(x) isnumeric(x));
    addRequired(ip, 'desc_prev',    @(x) isnumeric(x) && ismatrix(x));
    addRequired(ip, 'corners_prev', @(x) isnumeric(x));
    addParameter(ip, 'ncc_threshold', 0.7, @(x) isscalar(x));
    addParameter(ip, 'ratio_test',    0.8, @(x) isscalar(x) && x > 0 && x <= 1);
    addParameter(ip, 'max_distance',  50,  @(x) isscalar(x) && x > 0);
    addParameter(ip, 'mutual',        true, @islogical);
    parse(ip, desc_curr, corners_curr, desc_prev, corners_prev, varargin{:});

    ncc_thresh   = ip.Results.ncc_threshold;
    ratio_thresh = ip.Results.ratio_test;
    max_dist     = ip.Results.max_distance;
    do_mutual    = ip.Results.mutual;

    M = size(desc_curr, 1);
    N = size(desc_prev, 1);

    % Early exit
    if M == 0 || N == 0
        matches = zeros(0, 2);
        scores  = zeros(0, 1);
        return;
    end

    % ----------------------------------------------------------------
    % 1. Compute full NCC matrix via matrix multiply
    %    Since descriptors are zero-mean, unit-norm:
    %    NCC(i,j) = desc_curr(i,:) · desc_prev(j,:)'
    % ----------------------------------------------------------------
    NCC = desc_curr * desc_prev';   % [M x N]

    % ----------------------------------------------------------------
    % 2. Compute spatial distance matrix
    % ----------------------------------------------------------------
    dr = corners_curr(:,1) - corners_prev(:,1)';   % [M x N]
    dc = corners_curr(:,2) - corners_prev(:,2)';   % [M x N]
    dist_sq = dr.^2 + dc.^2;

    % ----------------------------------------------------------------
    % 3. Apply spatial gate
    % ----------------------------------------------------------------
    NCC(dist_sq > max_dist^2) = -Inf;

    % ----------------------------------------------------------------
    % 4. Forward matching: for each current corner, find best in prev
    % ----------------------------------------------------------------
    [sorted_ncc, sorted_idx] = sort(NCC, 2, 'descend');

    best_ncc     = sorted_ncc(:, 1);
    best_idx     = sorted_idx(:, 1);

    % Lowe's ratio test (need at least 2 candidates)
    if N >= 2
        second_ncc = sorted_ncc(:, 2);
        ratio_ok = (second_ncc <= 0) | (best_ncc ./ max(second_ncc, 1e-12) > 1 / ratio_thresh);
        % Explanation: ratio_test = best/second_best. We want
        % best/second < ratio_thresh to REJECT, so keep when
        % best/second >= 1/ratio_thresh... actually let me re-derive:
        % Lowe: reject if second_best / best > ratio_thresh
        % i.e. keep if second_best / best <= ratio_thresh
        ratio_ok = (second_ncc < 0) | ...
                   (second_ncc ./ max(best_ncc, 1e-12) <= ratio_thresh);
    else
        ratio_ok = true(M, 1);
    end

    % NCC threshold
    ncc_ok = best_ncc >= ncc_thresh;

    % Combined forward mask
    fwd_valid = ncc_ok & ratio_ok & isfinite(best_ncc);

    % ----------------------------------------------------------------
    % 5. Mutual matching (optional): check reverse direction
    % ----------------------------------------------------------------
    if do_mutual && any(fwd_valid)
        % For each prev corner, find best in curr
        [~, rev_best_idx] = max(NCC, [], 1);   % [1 x N]

        % A match (i -> j) is mutual if rev_best(j) == i
        for i = 1:M
            if fwd_valid(i)
                j = best_idx(i);
                if rev_best_idx(j) ~= i
                    fwd_valid(i) = false;
                end
            end
        end
    end

    % ----------------------------------------------------------------
    % 6. Assemble output
    % ----------------------------------------------------------------
    curr_indices = find(fwd_valid);
    prev_indices = best_idx(fwd_valid);
    match_scores = best_ncc(fwd_valid);

    % Sort by descending score
    [match_scores, sort_order] = sort(match_scores, 'descend');
    curr_indices = curr_indices(sort_order);
    prev_indices = prev_indices(sort_order);

    matches = [curr_indices(:), prev_indices(:)];
    scores  = match_scores(:);
end
