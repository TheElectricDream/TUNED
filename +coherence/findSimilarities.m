function [similarity_score, cv_map, regularity_map] = findSimilarities(sorted_x, sorted_y, iei_map, imgSz, radius)
% findSimilarities Compute local inter-event-interval regularity scores.
%
%   [similarity_score, cv_map, regularity_map] = findSimilarities( ...
%       sorted_x, sorted_y, iei_map, imgSz, radius)
%
%   Measures how temporally regular each event's spatial neighborhood is
%   by computing the local Coefficient of Variation (CV = sigma/mu) of the
%   persistent inter-event-interval (IEI) map. Real edges sweeping across
%   the sensor produce spatially coherent firing rates, yielding low CV.
%   Noise events fire at random rates relative to their neighbors, yielding
%   high CV. A pure-noise region where all pixels have similarly large IEIs
%   may still produce low CV, so the score also incorporates IEI magnitude
%   to disambiguate: only regions with both low CV AND reasonable IEI
%   magnitude receive high scores.
%
%   Inputs:
%     sorted_x  - Vector of x-coordinates (row indices) for events in the
%                 current frame, sorted by linear pixel index.
%     sorted_y  - Vector of y-coordinates (column indices) for events.
%     iei_map   - 2D matrix (imgSz) of persistent per-pixel inter-event
%                 intervals, e.g. maintained via an exponential moving
%                 average across frames. Pixels with no observed events
%                 should be 0.
%     imgSz     - Two-element vector [nrows, ncols].
%     radius    - (Optional) Scalar integer specifying the spatial
%                 neighborhood radius in pixels for the local CV
%                 computation. A radius of R produces a disk kernel of
%                 diameter (2R+1). Default: 4. Recommended range: 3-5.
%
%   Outputs:
%     similarity_score - N-by-1 vector of scores in [0, 1] for each input
%                        event. High values indicate temporally regular
%                        neighborhoods (likely real). Low values indicate
%                        irregular or inactive neighborhoods (likely noise).
%     cv_map           - 2D matrix (imgSz) of local coefficient of
%                        variation values. Useful for visualization and
%                        parameter tuning.
%     regularity_map   - 2D matrix (imgSz) of the final combined score
%                        before event lookup. Useful for visualization.
%
%   Algorithm:
%     1. Build an observation mask from pixels with nonzero IEI.
%     2. Use normalized convolution (disk kernel of given radius) to
%        compute local mean and local variance of IEI values, considering
%        only observed pixels.
%     3. Compute CV = local_std / local_mean at each pixel.
%     4. Convert CV to a regularity score in [0,1] via 1/(1 + CV).
%     5. Compute an IEI magnitude score that penalizes very large
%        intervals using a soft threshold at the median observed IEI.
%     6. Combine: score = regularity * magnitude_score.
%     7. Look up the combined score at each event location.


    % Validate input
    if nargin < 5
        radius = 4;
    end
    assert(isscalar(radius) && radius >= 1, ...
        'radius must be a positive integer >= 1.');
    assert(isequal(size(iei_map), imgSz), ...
        'iei_map dimensions must match imgSz.');

    % Cast to double for numerical stability in variance computation
    iei = double(iei_map);

    % Build observation mask and convolution kernel
    % Mask: 1 where we have IEI data, 0 elsewhere
    obs_mask = double(iei > 0);

    % Disk-shaped kernel for spatial neighborhood (isotropic weighting)
    kernel = fspecial('disk', radius);

    % Binary disk for counting observed neighbors (unweighted)
    count_kernel = double(fspecial('disk', radius) > 0);

    % Normalized convolution for local IEI statistics. Only observed events 
    % contribute to the local statistics.

    % Weighted sum of IEI values in the neighborhood
    weighted_sum = imfilter(iei .* obs_mask, kernel, 'replicate');

    % Sum of kernel weights over observed pixels
    weight_sum = imfilter(obs_mask, kernel, 'replicate');

    % Count of observed neighbors (unweighted, for minimum-count gate)
    neighbor_count = imfilter(obs_mask, count_kernel, 'replicate');

    % Local mean via normalized convolution
    local_mean = weighted_sum ./ max(weight_sum, eps);

    % Local variance via E[X^2] - E[X]^2
    weighted_sum_sq = imfilter((iei.^2) .* obs_mask, kernel, 'replicate');
    local_mean_sq = weighted_sum_sq ./ max(weight_sum, eps);
    local_var = max(local_mean_sq - local_mean.^2, 0);  % clamp rounding noise
    local_std = sqrt(local_var);

    % Coefficient of Variation
    cv_map = local_std ./ max(local_mean, eps);

    % Mark pixels with insufficient neighborhood data as invalid.
    % CV from fewer than 2 observed neighbors is not meaningful.
    insufficient_mask = (neighbor_count < 2) | (obs_mask == 0);
    cv_map(insufficient_mask) = Inf;  % maps to score = 0 below

    % Regularity score: low CV -> high score
    % Monotonic mapping from [0, inf) to (0, 1]:
    %   CV = 0   -> score = 1.0   (perfectly uniform neighborhood)
    %   CV = 1   -> score = 0.5   (std equals mean)
    %   CV -> inf -> score -> 0   (highly irregular)
    regularity_score = 1 ./ (1 + cv_map);

    % IEI magnitude score: penalize very large intervals. This addresses 
    % the failure mode where a pure-noise region has low CV
    % because all noise pixels fire at similarly slow rates. Such pixels
    % have large IEI values. We apply a soft penalty so that low-CV
    % regions only score well if their IEI is also reasonably small.
    %
    % The median observed IEI serves as a robust scale anchor:
    %   - At IEI = median,   score ~ 0.61
    %   - At IEI = 3*median, score ~ 0.22
    %   - At IEI << median,  score -> 1.0
    observed_iei = iei(obs_mask > 0);

    if ~isempty(observed_iei)
        iei_median = median(observed_iei);
        magnitude_score = exp(-iei ./ max(2 * iei_median, eps));
    else
        magnitude_score = zeros(imgSz);
    end

    % Inactive pixels get zero magnitude score
    magnitude_score(obs_mask == 0) = 0;

    % Combined regularity map
    regularity_map = regularity_score .* magnitude_score;

    % Look up per-event scores
    linear_idx = sub2ind(imgSz, sorted_x(:), sorted_y(:));
    similarity_score = regularity_map(linear_idx);

end