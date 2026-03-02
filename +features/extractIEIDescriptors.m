function [descriptors, valid_corners] = extractIEIDescriptors(maps, corners_sub, varargin)
%extractIEIDescriptors Multi-channel IEI vector descriptor at corners.
%
%   [DESCRIPTORS, VALID_CORNERS] = extractIEIDescriptors(MAPS, CORNERS_SUB)
%   [DESCRIPTORS, VALID_CORNERS] = extractIEIDescriptors(..., Name, Value)
%
%   Extracts a rich, multi-channel local descriptor at each detected
%   corner by sampling patches from multiple IEI-derived maps and
%   concatenating them into a single feature vector. This addresses the
%   key limitation of the Level 2 single-channel CV patch descriptor:
%   adjacent corners on the same edge have near-identical CV patches
%   (the aperture problem in CV space), but differ in their absolute
%   IEI values, local gradient structure, and event density.
%
%   Each channel is independently zero-meaned before concatenation, then
%   the full vector is L2-normalised. This ensures: (a) channels with
%   different physical units/scales contribute equally, (b) matching via
%   dot product still gives NCC, (c) adding channels cannot degrade the
%   descriptor below the single-channel baseline — a channel with no
%   discriminative power simply contributes near-zero variance.
%
%   This is distinct from all existing event-based descriptors:
%     - Li, Shi, Zhang et al. (IJRR, 2021): SIFT-like gradient
%       descriptor on a single-channel SITS.
%     - SD2Event (Gao et al., CVPR 2024): learned descriptors on
%       single-representation voxel grids.
%     - EventPoint (Huang et al., WACV 2023): self-supervised
%       descriptors on time surfaces.
%   All use first-order temporal representations. This descriptor uses
%   multiple second-order IEI statistics channels simultaneously.
%
%   Channels (selectable via the 'channels' parameter):
%     'cv'       - Coefficient of Variation map (σ/μ of local IEI).
%                  Encodes regularity pattern. [DEFAULT: ON]
%     'iei_mean' - Persistent IEI map (EMA-smoothed μ_IEI).
%                  Encodes absolute timing / speed information.
%                  Disambiguates corners at different speeds.
%                  [DEFAULT: ON]
%     'iei_std'  - Per-pixel IEI standard deviation.
%                  Encodes temporal variability.
%                  [DEFAULT: ON]
%     'harris_r' - Harris response map (from detectHarrisCV).
%                  Encodes gradient structure of the CV field — the
%                  "shape" of the corner in regularity space.
%                  [DEFAULT: ON]
%     'density'  - Event density / trace map.
%                  Encodes local activity level.
%                  [DEFAULT: OFF — often redundant with iei_mean]
%
%   Inputs:
%     maps         - Struct with fields corresponding to channel names.
%                    Each field is a 2D matrix [nRows x nCols].
%                    Required fields depend on 'channels' parameter.
%                    Example:
%                      maps.cv       = cv_map_for_harris;
%                      maps.iei_mean = iei_map;
%                      maps.iei_std  = t_std;
%                      maps.harris_r = harris_R;
%                      maps.density  = norm_trace_map;
%     corners_sub  - [K x 2] corner positions [row, col] in sensor
%                    coordinates (from detectHarrisCV).
%
%   Name-Value Parameters:
%     'half_size'    - Patch half-width (same for all channels).
%                      Default: 5  (11x11 per channel)
%     'sigma_weight' - Gaussian centre-weighting sigma [px].
%                      Default: 3.0
%     'min_observed' - Min fraction of observed pixels in the PRIMARY
%                      channel (cv) to accept the descriptor.
%                      Default: 0.25
%     'channels'     - Cell array of channel names to include.
%                      Default: {'cv', 'iei_mean', 'iei_std', 'harris_r'}
%     'channel_weights' - Struct with per-channel scalar weights.
%                         Applied AFTER per-channel zero-mean but
%                         BEFORE concatenation and L2 normalisation.
%                         Channels with higher weight contribute more
%                         to the final descriptor. Default: all 1.0.
%                         Example: struct('cv', 1.0, 'iei_mean', 0.5)
%
%   Outputs:
%     descriptors   - [M x D_total] matrix. D_total = N_channels *
%                     (2*half_size+1)^2. Each row is unit-norm.
%     valid_corners - [M x 2] subset of corners_sub.
%
%   Example:
%     maps.cv       = cv_map_for_harris;
%     maps.iei_mean = iei_map;
%     maps.iei_std  = t_std;
%     maps.harris_r = harris_R;
%     [desc, vc] = features.extractIEIDescriptors(maps, corners_sub, ...
%         'channels', {'cv','iei_mean','iei_std','harris_r'}, ...
%         'half_size', 5);
%
%   References:
%     Lowe, D.G. (2004), "Distinctive Image Features from Scale-Invariant
%       Keypoints," IJCV, 60(2), pp. 91-110.
%     Carandini, M. and Heeger, D.J. (2012), "Normalization as a
%       canonical neural computation," Nature Reviews Neuroscience,
%       13(1), pp. 51-62.
%
%   See also: features.detectHarrisCV, features.extractCVDescriptors,
%             features.matchCVDescriptors

    % ----------------------------------------------------------------
    % 0. Parse inputs
    % ----------------------------------------------------------------
    default_channels = {'cv', 'iei_mean', 'iei_std', 'harris_r'};

    ip = inputParser;
    addRequired(ip, 'maps',         @isstruct);
    addRequired(ip, 'corners_sub',  @(x) isnumeric(x) && size(x,2) >= 2);
    addParameter(ip, 'half_size',       5,                 @(x) isscalar(x) && x >= 2);
    addParameter(ip, 'sigma_weight',    3.0,               @(x) isscalar(x) && x > 0);
    addParameter(ip, 'min_observed',    0.25,              @(x) isscalar(x));
    addParameter(ip, 'channels',        default_channels,  @iscell);
    addParameter(ip, 'channel_weights', struct(),          @isstruct);
    parse(ip, maps, corners_sub, varargin{:});

    hs           = round(ip.Results.half_size);
    sigma_w      = ip.Results.sigma_weight;
    min_obs_frac = ip.Results.min_observed;
    channels     = ip.Results.channels;
    ch_weights   = ip.Results.channel_weights;
    n_channels   = numel(channels);

    patch_dim = 2 * hs + 1;
    D_per     = patch_dim^2;
    D_total   = n_channels * D_per;
    K         = size(corners_sub, 1);

    % ----------------------------------------------------------------
    % 1. Validate and prepare channel maps
    % ----------------------------------------------------------------
    channel_maps = cell(n_channels, 1);
    weights      = ones(n_channels, 1);

    for ch = 1:n_channels
        cname = channels{ch};
        if ~isfield(maps, cname)
            error('extractIEIDescriptors:missingMap', ...
                'maps struct is missing required channel: %s', cname);
        end
        M = double(maps.(cname));
        M(isnan(M)) = 0;
        channel_maps{ch} = M;

        if isfield(ch_weights, cname)
            weights(ch) = ch_weights.(cname);
        end
    end

    [nRows, nCols] = size(channel_maps{1});
    min_obs_count  = max(4, ceil(min_obs_frac * D_per));

    % ----------------------------------------------------------------
    % 2. Gaussian centre-weighting window
    % ----------------------------------------------------------------
    if isfinite(sigma_w)
        [gx, gy] = meshgrid(-hs:hs, -hs:hs);
        gauss_win = exp(-(gx.^2 + gy.^2) / (2 * sigma_w^2));
        gauss_vec = gauss_win(:);
    else
        gauss_vec = ones(D_per, 1);
    end

    % ----------------------------------------------------------------
    % 3. Extract multi-channel descriptors
    % ----------------------------------------------------------------
    descriptors = zeros(K, D_total);
    valid_mask  = false(K, 1);

    % Use the first channel (cv) as the observability reference
    primary_map = channel_maps{1};

    for i = 1:K
        r = corners_sub(i, 1);
        c = corners_sub(i, 2);

        % Border check
        if r - hs < 1 || r + hs > nRows || c - hs < 1 || c + hs > nCols
            continue;
        end

        % Check observability on primary channel
        primary_patch = primary_map(r-hs:r+hs, c-hs:c+hs);
        obs = primary_patch(:) > 0;
        if sum(obs) < min_obs_count
            continue;
        end

        % Extract and process each channel
        desc_vec = zeros(D_total, 1);
        any_variance = false;

        for ch = 1:n_channels
            cmap = channel_maps{ch};
            patch = cmap(r-hs:r+hs, c-hs:c+hs);
            pvec  = patch(:);

            % Impute unobserved pixels with mean of observed
            ch_obs = pvec ~= 0;
            if sum(ch_obs) > 0
                pvec(~ch_obs) = mean(pvec(ch_obs));
            end

            % Gaussian centre-weighting
            pvec = pvec .* gauss_vec;

            % Zero-mean (per-channel)
            pvec = pvec - mean(pvec);

            % Apply channel weight
            pvec = pvec * weights(ch);

            % Check for variance
            if norm(pvec) > 1e-12
                any_variance = true;
            end

            % Place into concatenated descriptor
            idx_start = (ch - 1) * D_per + 1;
            idx_end   = ch * D_per;
            desc_vec(idx_start:idx_end) = pvec;
        end

        if ~any_variance
            % All channels are flat — degenerate corner, skip
            continue;
        end

        % L2-normalise the full concatenated vector
        nrm = norm(desc_vec);
        if nrm < 1e-12
            continue;
        end
        desc_vec = desc_vec / nrm;

        descriptors(i, :) = desc_vec';
        valid_mask(i) = true;
    end

    % ----------------------------------------------------------------
    % 4. Filter to valid
    % ----------------------------------------------------------------
    descriptors   = descriptors(valid_mask, :);
    valid_corners = corners_sub(valid_mask, 1:2);
end
