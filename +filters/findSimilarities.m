function [similarity_score, cv_map, regularity_map] = findSimilarities(sorted_x, sorted_y, std_map, mean_map, imgSz)
% findSimilarities Compute per-pixel inter-event-interval regularity scores.
%
%   [similarity_score, cv_map, regularity_map] = findSimilarities( ...
%       sorted_x, sorted_y, std_map, mean_map, imgSz)
%
%   Measures how temporally regular each pixel's firing pattern is by
%   computing the Coefficient of Variation (CV = sigma/mu) directly from
%   the EMA-tracked per-pixel IEI statistics. Real signal pixels fire at
%   consistent intervals, yielding low CV. Noise pixels fire erratically,
%   yielding high CV. A pure-noise region where all pixels have similarly
%   large IEIs may still produce low CV, so the score also incorporates
%   IEI magnitude to disambiguate: only pixels with both low CV AND
%   reasonable mean IEI receive high scores.
%
%   Inputs:
%     sorted_x  - Vector of x-coordinates (row indices) for events in the
%                 current frame, sorted by linear pixel index.
%     sorted_y  - Vector of y-coordinates (column indices) for events.
%     std_map   - 2D matrix (imgSz) of EMA-tracked per-pixel IEI standard
%                 deviation. Pixels with no observed events should be 0.
%     mean_map  - 2D matrix (imgSz) of EMA-tracked per-pixel IEI mean.
%                 Pixels with no observed events should be 0.
%     imgSz     - Two-element vector [nrows, ncols].
%
%   Outputs:
%     similarity_score - N-by-1 vector of scores in [0, 1] for each input
%                        event. High values indicate temporally regular
%                        pixels (likely real). Low values indicate
%                        irregular or inactive pixels (likely noise).
%     cv_map           - 2D matrix (imgSz) of per-pixel coefficient of
%                        variation values. Useful for visualization,
%                        parameter tuning, and downstream feature
%                        detection (e.g. Harris on CV field).
%     regularity_map   - 2D matrix (imgSz) of the final combined score
%                        before event lookup. Useful for visualization.
%
%   Algorithm:
%     1. Build an observation mask from pixels with nonzero mean IEI.
%     2. Compute per-pixel CV = std_map / mean_map.
%     3. Convert CV to a regularity score in [0,1] via 1/(1 + CV).
%     4. Compute an IEI magnitude score that penalizes very large
%        intervals using a soft threshold at the median observed IEI.
%     5. Combine: score = regularity * magnitude_score.
%     6. Look up the combined score at each event location.
%
%   Notes:
%     - This function uses the pre-computed EMA-tracked per-pixel IEI
%       statistics (std_map, mean_map) rather than computing local
%       spatial statistics via convolution. Per-pixel temporal CV is the
%       principled metric for Rule 3: it directly measures whether a
%       pixel's own firing pattern is regular over time. Spatial context
%       is handled independently by Rule 1 (spatial density).
%     - The observation mask uses mean_map > 0 to identify pixels that
%       have received at least one event update.
%

    % ----------------------------------------------------------------
    % 1. Observation mask
    % ----------------------------------------------------------------
    obs_mask = mean_map > 0;

    % ----------------------------------------------------------------
    % 2. Per-pixel Coefficient of Variation
    % ----------------------------------------------------------------
    cv_map = std_map ./ max(mean_map, eps);

    % Unobserved pixels: set CV to Inf so they map to score = 0
    cv_map(~obs_mask) = Inf;

    % ----------------------------------------------------------------
    % 3. Regularity score: low CV -> high score
    %    Monotonic mapping from [0, inf) to (0, 1]:
    %      CV = 0   -> score = 1.0   (perfectly uniform firing)
    %      CV = 1   -> score = 0.5   (std equals mean)
    %      CV -> inf -> score -> 0   (highly irregular)
    % ----------------------------------------------------------------
    regularity_score = 1 ./ (1 + cv_map);

    % ----------------------------------------------------------------
    % 4. IEI magnitude score: penalize very large intervals
    %    Addresses the failure mode where a pure-noise region has
    %    low CV because all noise pixels fire at similarly slow rates.
    %    The median observed mean IEI serves as a robust scale anchor.
    % ----------------------------------------------------------------
    observed_mean_iei = mean_map(obs_mask);

    if ~isempty(observed_mean_iei)
        iei_median = median(observed_mean_iei);
        magnitude_score = exp(-mean_map ./ max(2 * iei_median, eps));
    else
        magnitude_score = zeros(imgSz);
    end

    % Inactive pixels get zero magnitude score
    magnitude_score(~obs_mask) = 0;

    % ----------------------------------------------------------------
    % 5. Combined regularity map
    % ----------------------------------------------------------------
    regularity_map = regularity_score .* magnitude_score;

    % ----------------------------------------------------------------
    % 6. Look up per-event scores
    % ----------------------------------------------------------------
    linear_idx = sub2ind(imgSz, sorted_x(:), sorted_y(:));
    similarity_score = regularity_map(linear_idx);

end