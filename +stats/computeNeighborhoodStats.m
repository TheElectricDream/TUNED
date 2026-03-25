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
%     All six outputs are computed without per-pixel loops by
%     exploiting the pre-sorted group structure of sorted_t:
%       1. min/max — direct indexing into group boundaries
%          (sorted_t is monotonic within each group).
%       2. mean — prefix-sum (cumsum) subtraction per group.
%       3. std  — prefix-sum of squared timestamps, then the
%          algebraic identity var(x) = E[x^2] - E[x]^2
%          with Bessel correction.
%       4. IEI mean — telescoping sum identity:
%          mean(diff(x)) = (x(end) - x(1)) / (n - 1).
%       5. IEI std — diff all of sorted_t in one pass, mask
%          out cross-group boundaries, then accumarray to
%          scatter sum(d^2) per group.
%
%   Notes:
%     - Input vectors must be pre-sorted by
%       stats.spreadEventsSpatially or equivalent grouping.
%     - Pixels with a single event get t_mean_diff = 0 and
%       t_std_diff = 0 (no interval observable).
%     - Pixels with two or fewer events get t_std_diff = 0
%       (need >= 2 IEIs for a sample standard deviation).
%     - The t_mean_diff output is the primary input to the
%       coherence IEI regularity rule and (via EMA smoothing)
%       the IEI-ATS accumulator tau mapping.
%     - The one-pass variance formula (E[x^2] - E[x]^2) can
%       suffer from catastrophic cancellation when absolute
%       values are large relative to the spread (Higham 2002,
%       Ch. 1). This is not an issue here because timestamps
%       are relative within each frame window, keeping
%       magnitudes moderate. A max(..., 0) guard clamps any
%       negative values caused by floating-point noise.
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

    K = numel(unique_idx);
    if K == 0
        return;
    end

    % ----------------------------------------------------------------
    % 1. Group counts
    % ----------------------------------------------------------------
    counts = group_ends - pos + 1;         % [K x 1] events per pixel

    % ----------------------------------------------------------------
    % 2. Min / max — direct indexing (sorted within group)
    % ----------------------------------------------------------------
    mn = sorted_t(pos);
    mx = sorted_t(group_ends);

    % ----------------------------------------------------------------
    % 3. Mean — prefix-sum subtraction
    % ----------------------------------------------------------------
    % Prepend zero so that cs(end+1) - cs(start) gives the group sum
    % without special-casing the first group.
    cs      = [0; cumsum(sorted_t)];
    grp_sum = cs(group_ends + 1) - cs(pos);
    mu      = grp_sum ./ counts;

    % ----------------------------------------------------------------
    % 4. Std — prefix-sum of squares + algebraic variance
    % ----------------------------------------------------------------
    %   var(x) = [ sum(x^2) - n * mu^2 ] / (n - 1)
    cs2      = [0; cumsum(sorted_t .^ 2)];
    grp_sum2 = cs2(group_ends + 1) - cs2(pos);

    denom    = max(counts - 1, 1);         % protect single-event pixels
    variance = (grp_sum2 - counts .* mu.^2) ./ denom;
    sd       = sqrt(max(variance, 0));     % clamp float noise
    sd(counts == 1) = 0;

    % ----------------------------------------------------------------
    % 5. IEI mean — telescoping sum identity
    % ----------------------------------------------------------------
    %   mean(diff(x)) = (x(end) - x(1)) / (numel(x) - 1)
    diff_counts = counts - 1;
    mu_d = (mx - mn) ./ max(diff_counts, 1);
    mu_d(counts <= 1) = 0;

    % ----------------------------------------------------------------
    % 6. IEI std — vectorized diff with cross-group masking
    % ----------------------------------------------------------------
    M     = numel(sorted_t);
    d_all = diff(sorted_t);               % [M-1 x 1]

    % 6a. Build per-element group labels via cumsum delta trick:
    %     place a 1 at every group start, cumsum maps each element
    %     to its group index.
    delta     = zeros(M, 1);
    delta(pos) = 1;
    g         = cumsum(delta);            % g(i) = group index of sorted_t(i)

    % 6b. Identify within-group diffs (exclude boundary diffs where
    %     consecutive elements belong to different groups).
    g_diff = g(1:end-1);                  % group label for each diff
    within = true(M - 1, 1);
    if K > 1
        within(group_ends(1:end-1)) = false;
    end

    % 6c. Scatter sum(d^2) per group via accumarray
    d_w    = d_all(within);
    gw     = g_diff(within);
    sum_d2 = accumarray(gw, d_w .^ 2, [K, 1]);

    %   var(d) = [ sum(d^2)/n - mu_d^2 ] * n/(n-1)    (Bessel)
    n_iei = max(diff_counts, 1);
    var_d = (sum_d2 ./ n_iei - mu_d .^ 2) .* n_iei ./ max(n_iei - 1, 1);
    sd_d  = sqrt(max(var_d, 0));
    sd_d(counts <= 2) = 0;               % need >= 2 IEIs for std

    % ----------------------------------------------------------------
    % 7. Scatter into image maps
    % ----------------------------------------------------------------
    t_mean(unique_idx)      = mu;
    t_std(unique_idx)       = sd;
    t_max(unique_idx)       = mx;
    t_min(unique_idx)       = mn;
    t_mean_diff(unique_idx) = mu_d;
    t_std_diff(unique_idx)  = sd_d;

end