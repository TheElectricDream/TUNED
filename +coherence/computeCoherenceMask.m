function [norm_trace_map, norm_similarity_map, ...
    norm_persist_map, filtered_coherence_map] = ...
    computeCoherenceMask(sorted_x, sorted_y, sorted_t, imgSz, ...
    t_interval, unique_idx, pos, group_ends, coh_params, ...
    frameIndex, norm_trace_map_prev, iei_map)
% COMPUTECOHERENCEMASK  Three-rule coherence filtering pipeline.
%
%   [NORM_TRACE_MAP, NORM_SIMILARITY_MAP, NORM_PERSIST_MAP,
%   FILTERED_COHERENCE_MAP] = COMPUTECOHERENCEMASK(...) combines
%   three independent coherence rules into a scalar per-pixel score
%   that separates real edge events from noise before temporal
%   surface accumulation.
%
%   Inputs:
%     sorted_x, sorted_y  - [N x 1] Pixel coordinates (row, col),
%                           sorted by linear pixel index.
%     sorted_t             - [N x 1] Timestamps [s], sorted to match.
%     imgSz                - [1 x 2] Image dimensions [nRows, nCols].
%     t_interval           - Scalar frame interval [s].
%     unique_idx           - [K x 1] Linear indices of active pixels.
%     pos                  - [K x 1] Start of each pixel group in the
%                           sorted vectors.
%     group_ends           - [K x 1] End of each pixel group.
%     coh_params           - Struct with fields:
%       .r_s                   - Spatial search radius (normalized)
%       .trace_threshold       - Min neighbour distance sum for Rule 1
%       .persistence_threshold - Max cross-frame distance for Rule 2
%       .similarity_threshold  - Max regularity score for Rule 3
%     frameIndex           - Current frame number (1-indexed). Rule 2
%                           is skipped on the first frame.
%     norm_trace_map_prev  - [imgSz] Previous frame's trace map for
%                           persistence assessment (Rule 2).
%     iei_map              - [imgSz] Per-pixel IEI map (current-frame
%                           mean diff or EMA-smoothed) for Rule 3.
%
%   Outputs:
%     norm_trace_map       - [imgSz] Rule 1: spatial density score.
%     norm_similarity_map  - [imgSz] Rule 3: IEI regularity score.
%     norm_persist_map     - [imgSz] Rule 2: temporal persistence.
%     filtered_coherence_map - [imgSz] Combined score (sum of rules).
%
%   Algorithm:
%     1. Rule 1 — Spatial density: KD-tree radius search over
%        normalized (x, y, t) space. Sum of distances within r_s
%        forms a density proxy. Log-normalized to compress the
%        heavy-tailed distribution.
%     2. Rule 3 — IEI regularity: local coefficient of variation
%        (CV = sigma/mu) of the IEI map via normalized convolution.
%        Scores above similarity_threshold are rejected.
%     3. Rule 2 — Temporal persistence: KNN search between the
%        current and previous trace maps in normalized (row, col,
%        value) space. Pixels without a close predecessor are
%        rejected. Skipped on frame 1.
%     4. Combine: sum of the three normalized rule maps.
%
%   Notes:
%     - The three rules are independent and additive. Each produces
%       a [0, 1] normalized map. The combined map is thresholded
%       downstream in main.m (not inside this function).
%     - Coordinates: x = row, y = col, sub2ind(imgSz, x, y).
%
%   See also: coherence.findSpatialNeighbours,
%             coherence.findSimilarities,
%             coherence.findPersistenceVectorized

    % ----------------------------------------------------------------
    % 0. Parse parameters
    % ----------------------------------------------------------------
    r_s                   = coh_params.r_s;
    trace_threshold       = coh_params.trace_threshold;
    persistence_threshold = coh_params.persistence_threshold;
    similarity_threshold  = coh_params.similarity_threshold;

    % ----------------------------------------------------------------
    % 1. Rule 1 — Spatial density (trace map)
    % ----------------------------------------------------------------
    sum_exp_dist_map = zeros(imgSz);

    % KD-tree radius search in normalized (x, y, t) space
    [~, distances_db] = coherence.findSpatialNeighbours(...
        sorted_x, sorted_y, sorted_t, r_s, imgSz, t_interval);

    % Sum distances for each event
    sum_exp_event = cellfun(@sum, distances_db);

    % Aggregate to pixel map (max over events at same pixel)
    for k = 1:length(unique_idx)
        val_chunk_exp = sum_exp_event(pos(k):group_ends(k));
        idx = unique_idx(k);
        sum_exp_dist_map(idx) = max(val_chunk_exp);
    end

    % Threshold: reject spatially isolated events
    sum_exp_dist_map(sum_exp_dist_map <= trace_threshold) = 0;

    % Log-normalize to compress the heavy-tailed distribution
    log_trace_map = log1p(sum_exp_dist_map');
    norm_trace_map = log_trace_map' ./ max(log_trace_map(:));

    % ----------------------------------------------------------------
    % 2. Rule 3 — IEI regularity (similarity map)
    % ----------------------------------------------------------------
    [~, ~, norm_similarity_map] = coherence.findSimilarities(...
        sorted_x, sorted_y, iei_map, imgSz, 10);

    % Reject regions with regularity score above threshold
    norm_similarity_map(norm_similarity_map > similarity_threshold) ...
        = nan;
    norm_similarity_map(isnan(norm_similarity_map)) = 0;

    % ----------------------------------------------------------------
    % 3. Rule 2 — Temporal persistence
    % ----------------------------------------------------------------
    if frameIndex == 1
        % No previous frame available — use trace map as proxy
        persist_map = norm_trace_map;
    else
        persist_map = zeros(size(norm_trace_map));

        % KNN search between current and previous trace maps
        [~, ~, minDists, validIdx] = ...
            coherence.findPersistenceVectorized(...
            norm_trace_map, norm_trace_map_prev, imgSz);

        if ~isempty(validIdx)
            persist_map(validIdx) = minDists;
        end

        % Reject pixels without a close predecessor
        persist_map(persist_map > persistence_threshold) = nan;
    end

    % Log-normalize the persistence map
    log_persist_map = log1p(persist_map);
    norm_persist_map = log_persist_map ./ max(log_persist_map(:));

    % ----------------------------------------------------------------
    % 4. Combine rule maps
    % ----------------------------------------------------------------
    filtered_coherence_map = norm_trace_map ...
        + norm_persist_map + norm_similarity_map;

end