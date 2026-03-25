function R = computeClarkEvans(map, varargin)
% COMPUTECLARKEVANS  Clark-Evans nearest-neighbour aggregation index.
%
%   R = COMPUTECLARKEVANS(MAP) extracts active pixel coordinates from
%   MAP (map > 0) and computes the Clark-Evans R statistic — the ratio
%   of the observed mean nearest-neighbour distance to the distance
%   expected under complete spatial randomness (CSR).
%
%   R = COMPUTECLARKEVANS(MAP, 'Name', Value) accepts optional
%   name-value arguments:
%     'ImageSize' - ([H W], derived from MAP) Area of the observation
%                   window in pixels. Defaults to size(map), so
%                   override only when MAP is a cropped ROI of a
%                   larger sensor.
%     'MaxPoints' - (0) Maximum number of active pixels to include
%                   in the KD-tree query. When N > MaxPoints, a
%                   uniform random subsample of size MaxPoints is
%                   drawn. The expected NND under CSR is computed
%                   using the subsampled count, so the R ratio
%                   remains consistent. Set to 0 (default) to
%                   disable subsampling and use all active pixels.
%
%   Inputs:
%     map - [H x W] 2D value map. Active pixels are those with
%           map > 0.
%
%   Outputs:
%     R - Scalar Clark-Evans index.
%           R ≈ 1 : randomly distributed (CSR)
%           R < 1 : clustered
%           R > 1 : over-dispersed / regular
%
%   Notes:
%     - Nearest-neighbour search uses KNNSEARCH (k = 2,
%       self-excluded) rather than an O(N^2) explicit loop, giving
%       O(N log N) scaling.
%     - When MaxPoints is active, the statistic is computed over the
%       subsample. Both r_obs and r_exp use the subsampled count N_q,
%       so R is an unbiased estimator of the full-set R at large N.
%       The observation window area A remains the full image area,
%       since the subsample is drawn from the full spatial extent.
%     - Returns NaN and emits a warning when fewer than 2 active
%       pixels are present, as the statistic is undefined.
%
%   Reference:
%     P. J. Clark and F. C. Evans, "Distance to Nearest Neighbor as
%     a Measure of Spatial Relationships in Populations," Ecology,
%     vol. 35, no. 4, pp. 445-453, Oct. 1954.
%     DOI: 10.2307/1931034
%
%   See also: plot.showSpreadDistribution, stats.findSpatialNeighbours

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    [H, W] = size(map);
    p = inputParser;
    addParameter(p, 'ImageSize', [H, W], ...
        @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));
    addParameter(p, 'MaxPoints', 0, ...
        @(x) isscalar(x) && isnumeric(x) && x >= 0);
    parse(p, varargin{:});
    img_sz    = p.Results.ImageSize;
    max_pts   = round(p.Results.MaxPoints);

    % ----------------------------------------------------------------
    % 1. Extract active pixel coordinates [col, row]
    % ----------------------------------------------------------------
    [rows, cols] = find(map > 0);
    N_full = numel(rows);

    if N_full < 2
        warning('computeClarkEvans:insufficientPoints', ...
            'Fewer than 2 active pixels — R statistic is undefined.');
        R = NaN;
        return;
    end

    % ----------------------------------------------------------------
    % 2. (Optional) Subsample to cap KD-tree cost
    % ----------------------------------------------------------------
    if max_pts > 0 && N_full > max_pts
        idx_sub = randperm(N_full, max_pts);
        pts = [cols(idx_sub), rows(idx_sub)];
        N_q = max_pts;
    else
        pts = [cols, rows];
        N_q = N_full;
    end

    % ----------------------------------------------------------------
    % 3. Nearest-neighbour distances  (O(N_q log N_q) via kd-tree)
    % ----------------------------------------------------------------
    % k = 2: first neighbour is the point itself; second is the
    % nearest other point.
    [~, D_nn] = knnsearch(pts, pts, 'K', 2);
    r_obs = mean(D_nn(:, 2));

    % ----------------------------------------------------------------
    % 4. Clark-Evans statistic
    % ----------------------------------------------------------------
    A     = img_sz(1) * img_sz(2);     % observation window area [px^2]
    r_exp = 0.5 / sqrt(N_q / A);      % expected NND under CSR
    R     = r_obs / r_exp;

end