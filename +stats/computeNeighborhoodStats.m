function [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
    computeNeighborhoodStats(sorted_t, unique_idx, pos, ...
    group_ends, imgSz)
% COMPUTENEIGHBORHOODSTATS  Per-pixel IEI statistics from grouped events.
%
%   [T_MEAN, T_STD, T_MAX, T_MIN, T_MEAN_DIFF, T_STD_DIFF] =
%   COMPUTENEIGHBORHOODSTATS(SORTED_T, UNIQUE_IDX, POS, GROUP_ENDS,
%   IMGSZ) computes timestamp statistics for each pixel that received
%   events in the current time window. The inter-event interval (IEI)
%   statistics (t_mean_diff, t_std_diff) are derived from consecutive
%   timestamp differences within each pixel's event sequence.
%
%   Inputs:
%     sorted_t   - [M x 1] Timestamps sorted by linear pixel index,
%                  then by time within each pixel group.
%     unique_idx - [K x 1] Linear pixel indices for each active pixel.
%     pos        - [K x 1] Start position of each pixel's group in
%                  sorted_t.
%     group_ends - [K x 1] End position of each pixel's group in
%                  sorted_t.
%     imgSz      - [1 x 2] Image dimensions [nRows, nCols].
%
%   Outputs:
%     t_mean      - [imgSz] Mean timestamp per pixel.
%     t_std       - [imgSz] Standard deviation of timestamps.
%     t_max       - [imgSz] Maximum timestamp per pixel.
%     t_min       - [imgSz] Minimum timestamp per pixel.
%     t_mean_diff - [imgSz] Mean inter-event interval (IEI) per pixel.
%                   Zero if the pixel received only one event.
%     t_std_diff  - [imgSz] Standard deviation of IEI per pixel.
%
%   Algorithm:
%     For each active pixel:
%       1. Extract the timestamp subsequence.
%       2. Compute mean, std, max, min of raw timestamps.
%       3. If more than one event, compute diff() and derive
%          mean and std of the inter-event intervals.
%
%   Notes:
%     - Input vectors must be pre-sorted by
%       stats.spreadEventsSpatially or equivalent grouping.
%     - Pixels with a single event get t_mean_diff = 0 and
%       t_std_diff = 0 (no interval observable).
%     - The t_mean_diff output is the primary input to the
%       coherence IEI regularity rule and (via EMA smoothing)
%       the IEI-ATS accumulator tau mapping.
%
%   See also: stats.spreadEventsSpatially,
%             coherence.computeCoherenceMask

    % ----------------------------------------------------------------
    % 0. Initialize output maps
    % ----------------------------------------------------------------
    t_mean      = zeros(imgSz);
    t_std       = zeros(imgSz);
    t_max       = zeros(imgSz);
    t_min       = zeros(imgSz);
    t_mean_diff = zeros(imgSz);
    t_std_diff  = zeros(imgSz);

    % ----------------------------------------------------------------
    % 1. Compute per-pixel statistics
    % ----------------------------------------------------------------
    for k = 1:length(unique_idx)
        val_chunk_t = sorted_t(pos(k):group_ends(k));
        idx = unique_idx(k);

        % Raw timestamp statistics
        t_mean(idx) = mean(val_chunk_t);
        t_max(idx)  = max(val_chunk_t);
        t_min(idx)  = min(val_chunk_t);
        t_std(idx)  = std(val_chunk_t);

        % Inter-event interval statistics
        if numel(val_chunk_t) > 1
            d = diff(val_chunk_t);
            t_mean_diff(idx) = mean(d);
            t_std_diff(idx)  = std(d);
        else
            t_mean_diff(idx) = 0;
            t_std_diff(idx)  = 0;
        end
    end

end