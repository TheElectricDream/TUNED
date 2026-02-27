function [descriptors, valid_corners] = extractCVDescriptors(cv_map, corners_sub, varargin)
%extractCVDescriptors Extract local CV-patch descriptors at corner locations.
%
%   [DESCRIPTORS, VALID_CORNERS] = extractCVDescriptors(CV_MAP, CORNERS_SUB)
%   [DESCRIPTORS, VALID_CORNERS] = extractCVDescriptors(..., Name, Value)
%
%   For each detected corner, extracts a local patch from the Coefficient
%   of Variation (CV) map and converts it into a normalised descriptor
%   vector. The descriptor encodes the *local pattern of temporal
%   regularity* around a corner — a second-order IEI statistic that is
%   inherently less sensitive to absolute event rate (and therefore to
%   feature speed) than timestamp-based descriptors.
%
%   This is distinct from all existing event-based descriptors:
%     - Li, Shi, Zhang et al. (IJRR, 2021) build SIFT-like gradient
%       descriptors on the Speed-Invariant Time Surface (SITS).
%     - SD2Event (Gao et al., CVPR 2024) learns descriptors on voxel
%       grids.
%     - EventPoint (Huang et al., WACV 2023) uses self-supervised
%       descriptors on time surfaces.
%   All operate on first-order temporal representations. This descriptor
%   operates on second-order IEI statistics (σ/μ).
%
%   Algorithm:
%     1. Extract a (2*half_size+1)^2 patch from cv_map at each corner.
%     2. Impute unobserved pixels (NaN/0) with patch mean of observed
%        pixels (neutral imputation — no bias toward any CV value).
%     3. Apply Gaussian centre-weighting (reduces sensitivity to
%        localisation error, following the SIFT philosophy).
%     4. Zero-mean, L2-normalise → unit-norm descriptor.
%        Matching via dot product gives NCC directly.
%     5. Corners too close to the border are discarded.
%
%   Inputs:
%     cv_map       - 2D matrix of local CV values (output of
%                    coherence.findSimilarities).
%     corners_sub  - [K x 2] array of [row, col] in sensor coordinates
%                    (from detectHarrisCV).
%
%   Name-Value Parameters:
%     'half_size'    - Half-width of descriptor patch.
%                      Full patch: (2*half_size+1) x (2*half_size+1).
%                      Default: 5  → 11x11 = 121-dim descriptor.
%     'sigma_weight' - Gaussian centre-weighting sigma [px].
%                      Set to Inf for uniform weighting.
%                      Default: 3.0
%     'min_observed' - Minimum fraction of observed (nonzero) pixels in
%                      the patch to accept the descriptor. [0, 1].
%                      Default: 0.25
%
%   Outputs:
%     descriptors   - [M x D] matrix. M ≤ K valid descriptors, each
%                     D = (2*half_size+1)^2. Rows are unit-norm,
%                     zero-mean vectors. Dot product between two
%                     descriptors = NCC score in [-1, 1].
%     valid_corners - [M x 2] subset of corners_sub that produced
%                     valid descriptors.
%
%   Example:
%     [desc, vc] = features.extractCVDescriptors(cv_map, corners_sub, ...
%         'half_size', 5, 'sigma_weight', 3.0);
%     % desc(i,:) · desc(j,:)' gives NCC between corners i and j.
%
%   References:
%     Lowe, D.G. (2004), "Distinctive Image Features from Scale-Invariant
%       Keypoints," IJCV, 60(2), pp. 91–110.
%
%   See also: features.detectHarrisCV, features.matchCVDescriptors

    % ----------------------------------------------------------------
    % 0. Parse inputs
    % ----------------------------------------------------------------
    ip = inputParser;
    addRequired(ip, 'cv_map',      @(x) isnumeric(x) && ismatrix(x));
    addRequired(ip, 'corners_sub', @(x) isnumeric(x) && size(x,2) >= 2);
    addParameter(ip, 'half_size',    5,    @(x) isscalar(x) && x >= 2);
    addParameter(ip, 'sigma_weight', 3.0,  @(x) isscalar(x) && x > 0);
    addParameter(ip, 'min_observed', 0.25, @(x) isscalar(x) && x >= 0 && x <= 1);
    parse(ip, cv_map, corners_sub, varargin{:});

    hs          = round(ip.Results.half_size);
    sigma_w     = ip.Results.sigma_weight;
    min_obs_frac = ip.Results.min_observed;

    [nRows, nCols] = size(cv_map);
    patch_dim = 2 * hs + 1;
    D = patch_dim^2;
    K = size(corners_sub, 1);
    min_obs_count = max(4, ceil(min_obs_frac * D));

    % ----------------------------------------------------------------
    % 1. Gaussian centre-weighting window
    % ----------------------------------------------------------------
    if isfinite(sigma_w)
        [gx, gy] = meshgrid(-hs:hs, -hs:hs);
        gauss_win = exp(-(gx.^2 + gy.^2) / (2 * sigma_w^2));
        gauss_vec = gauss_win(:);
    else
        gauss_vec = ones(D, 1);
    end

    % ----------------------------------------------------------------
    % 2. Sanitise CV map
    % ----------------------------------------------------------------
    S = double(cv_map);
    S(isnan(S)) = 0;

    % ----------------------------------------------------------------
    % 3. Extract, weight, normalise
    % ----------------------------------------------------------------
    descriptors   = zeros(K, D);
    valid_mask    = false(K, 1);

    for i = 1:K
        r = corners_sub(i, 1);
        c = corners_sub(i, 2);

        % Border check
        if r - hs < 1 || r + hs > nRows || c - hs < 1 || c + hs > nCols
            continue;
        end

        % Extract raw patch
        patch = S(r-hs:r+hs, c-hs:c+hs);
        pvec  = patch(:);

        % Count observed pixels
        obs = pvec > 0;
        n_obs = sum(obs);
        if n_obs < min_obs_count
            continue;
        end

        % Impute unobserved with mean of observed (neutral fill)
        if n_obs < D
            pvec(~obs) = mean(pvec(obs));
        end

        % Gaussian centre-weighting
        pvec = pvec .* gauss_vec;

        % Zero-mean
        pvec = pvec - mean(pvec);

        % L2 normalise
        nrm = norm(pvec);
        if nrm < 1e-12
            % Zero-variance patch — degenerate, skip
            continue;
        end
        pvec = pvec / nrm;

        descriptors(i, :) = pvec';
        valid_mask(i) = true;
    end

    % ----------------------------------------------------------------
    % 4. Filter to valid
    % ----------------------------------------------------------------
    descriptors   = descriptors(valid_mask, :);
    valid_corners = corners_sub(valid_mask, 1:2);
end
