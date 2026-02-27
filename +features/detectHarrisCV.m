function [corners, R, corners_sub] = detectHarrisCV(cv_map, varargin)
%detectHarrisCV Detect corners in a Coefficient of Variation map.
%
%   [CORNERS, R, CORNERS_SUB] = detectHarrisCV(CV_MAP)
%   [CORNERS, R, CORNERS_SUB] = detectHarrisCV(CV_MAP, Name, Value, ...)
%
%   Applies the Harris corner response to a 2D Coefficient of Variation
%   (CV) map derived from inter-event-interval (IEI) statistics. Corners
%   in the CV field correspond to locations where the *temporal regularity
%   landscape* changes in two spatial directions simultaneously — i.e.,
%   where the local pattern of event firing regularity is distinctive.
%
%   This is distinct from all existing event-based Harris variants
%   (eHarris, FA-Harris, luvHarris, CHEC, SE-Harris, DTFS-eHarris),
%   which operate on first-order timestamp surfaces (SAE or decay state).
%   The CV map encodes second-order temporal statistics (σ/μ of IEI),
%   making the detector inherently more invariant to absolute event rate
%   and, by extension, to feature speed.
%
%   Algorithm:
%     1. Compute spatial gradients of the CV map via Sobel operators.
%     2. Form the structure tensor components (Ix^2, Iy^2, Ix*Iy).
%     3. Smooth the structure tensor with a Gaussian kernel of width
%        sigma_smooth (integrates gradient information over a local
%        neighbourhood — the "auto-correlation" window).
%     4. Compute Harris response:  R = det(M) - k * trace(M)^2
%        where M = [sum(Ix^2)  sum(Ix*Iy);  sum(Ix*Iy)  sum(Iy^2)]
%     5. Threshold R and apply non-maximum suppression (NMS).
%
%   Inputs:
%     cv_map     - 2D matrix (same size as sensor) of local CV values.
%                  Output of coherence.findSimilarities (the cv_map
%                  output) or any map where low values indicate regular
%                  regions and high values indicate irregular ones.
%                  NaN and zero entries are treated as unobserved.
%
%   Name-Value Parameters:
%     'k'            - Harris sensitivity parameter. Typical range
%                      [0.04, 0.15]. Larger values suppress edges more
%                      aggressively, retaining only strong corners.
%                      Default: 0.06
%     'sigma_smooth' - Gaussian smoothing sigma for the structure tensor
%                      (in pixels). Controls the spatial integration
%                      scale. Larger values detect coarser corners.
%                      Default: 1.5
%     'threshold'    - Fraction of max(R) below which responses are
%                      suppressed. Range [0, 1]. Default: 0.01
%     'nms_radius'   - Radius (pixels) for non-maximum suppression.
%                      Default: 5
%     'max_corners'  - Maximum number of corners to return (strongest
%                      first). Set to Inf for no limit. Default: 200
%     'border'       - Border exclusion zone in pixels. Corners within
%                      this margin of the image edge are discarded.
%                      Default: 5
%
%   Outputs:
%     corners     - [K x 3] array of [row, col, response] for each
%                   detected corner, sorted by descending response.
%     R           - 2D Harris response map (same size as cv_map).
%                   Useful for visualization and parameter tuning.
%     corners_sub - [K x 2] array of [row, col] only (convenience
%                   output for overlay plotting).
%
%   Example:
%     % After computing the CV map from findSimilarities:
%     [~, cv_map, ~] = coherence.findSimilarities(sx, sy, iei_map, imgSz);
%     [corners, R] = features.detectHarrisCV(cv_map, ...
%         'k', 0.06, 'nms_radius', 5, 'threshold', 0.01);
%     % Overlay on existing figure:
%     hold on;
%     plot(corners(:,2), corners(:,1), 'r+', 'MarkerSize', 8);
%
%   Motivation:
%     Existing event-based Harris detectors compute the structure tensor
%     of a Surface of Active Events (timestamps). The CV map instead
%     encodes the *regularity* of event arrival — a second-order
%     statistic that is inherently less sensitive to absolute speed.
%     Corners in the CV field are locations where the temporal regularity
%     pattern is locally distinctive, making them candidates for
%     speed-invariant feature detection in mixed-dynamics scenarios
%     (e.g., tumbling spacecraft with fast-moving edges and slow
%     background stars).
%
%   References:
%     Harris, C. and Stephens, M. (1988), "A Combined Corner and Edge
%       Detector," Proc. 4th Alvey Vision Conference, pp. 147–151.
%     Abdi, H. (2010), "Coefficient of Variation," in Encyclopedia of
%       Research Design, SAGE Publications.
%     Knutsson, H. and Westin, C.-F. (1993), "Normalized and
%       Differential Convolution," Proc. IEEE CVPR, pp. 515–523.
%
%   See also: coherence.findSimilarities, features.detectArcStarCorners

    % ----------------------------------------------------------------
    % 0. Parse inputs
    % ----------------------------------------------------------------
    p = inputParser;
    addRequired(p, 'cv_map', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'k',            0.06,  @(x) isscalar(x) && x > 0);
    addParameter(p, 'sigma_smooth', 1.5,   @(x) isscalar(x) && x > 0);
    addParameter(p, 'threshold',    0.01,  @(x) isscalar(x) && x >= 0 && x <= 1);
    addParameter(p, 'nms_radius',   5,     @(x) isscalar(x) && x >= 1);
    addParameter(p, 'max_corners',  200,   @(x) isscalar(x) && x >= 1);
    addParameter(p, 'border',       5,     @(x) isscalar(x) && x >= 0);
    parse(p, cv_map, varargin{:});

    k            = p.Results.k;
    sigma_smooth = p.Results.sigma_smooth;
    thresh_frac  = p.Results.threshold;
    nms_r        = round(p.Results.nms_radius);
    max_corners  = p.Results.max_corners;
    border       = round(p.Results.border);

    [nRows, nCols] = size(cv_map);

    % ----------------------------------------------------------------
    % 1. Sanitise input: replace NaN/zero with 0 (unobserved)
    % ----------------------------------------------------------------
    S = double(cv_map);
    S(isnan(S)) = 0;

    % ----------------------------------------------------------------
    % 2. Spatial gradients via Sobel operators
    % ----------------------------------------------------------------
    % Standard 3x3 Sobel kernels
    Kx = [-1 0 1; -2 0 2; -1 0 1] / 8;   % normalised
    Ky = [-1 -2 -1;  0  0  0;  1  2  1] / 8;

    Ix = imfilter(S, Kx, 'replicate');
    Iy = imfilter(S, Ky, 'replicate');

    % ----------------------------------------------------------------
    % 3. Structure tensor components, smoothed with Gaussian
    % ----------------------------------------------------------------
    % The Gaussian integration window is what makes Harris robust —
    % it aggregates gradient information over a spatial neighbourhood.
    Ix2  = imgaussfilt(Ix .* Ix, sigma_smooth);
    Iy2  = imgaussfilt(Iy .* Iy, sigma_smooth);
    IxIy = imgaussfilt(Ix .* Iy, sigma_smooth);

    % ----------------------------------------------------------------
    % 4. Harris response:  R = det(M) - k * trace(M)^2
    % ----------------------------------------------------------------
    detM   = Ix2 .* Iy2 - IxIy .* IxIy;
    traceM = Ix2 + Iy2;
    R      = detM - k * (traceM .^ 2);

    % Zero out unobserved regions (where CV was 0/NaN)
    R(S == 0) = 0;

    % ----------------------------------------------------------------
    % 5. Threshold
    % ----------------------------------------------------------------
    R_max = max(R(:));
    if R_max <= 0
        corners     = zeros(0, 3);
        corners_sub = zeros(0, 2);
        return;
    end
    R_thresh = R;
    R_thresh(R < thresh_frac * R_max) = 0;

    % ----------------------------------------------------------------
    % 6. Border exclusion
    % ----------------------------------------------------------------
    if border > 0
        R_thresh(1:border, :)           = 0;
        R_thresh(end-border+1:end, :)   = 0;
        R_thresh(:, 1:border)           = 0;
        R_thresh(:, end-border+1:end)   = 0;
    end

    R_thresh(isnan(R_thresh))=0;

    % ----------------------------------------------------------------
    % 7. Non-maximum suppression via morphological dilation
    % ----------------------------------------------------------------
    % A pixel is a local maximum if it equals the max in its
    % neighbourhood after dilation.
    se       = strel('disk', nms_r, 0);
    R_dilate = imdilate(R_thresh, se);
    nms_mask = (R_thresh == R_dilate) & (R_thresh > 0);

    % ----------------------------------------------------------------
    % 8. Extract and sort corners
    % ----------------------------------------------------------------
    [rows, cols] = find(nms_mask);
    if isempty(rows)
        corners     = zeros(0, 3);
        corners_sub = zeros(0, 2);
        return;
    end

    idx       = sub2ind([nRows, nCols], rows, cols);
    responses = R(idx);

    % Sort by descending response strength
    [responses, sort_idx] = sort(responses, 'descend');
    rows = rows(sort_idx);
    cols = cols(sort_idx);

    % Limit to max_corners
    n_keep = min(numel(rows), max_corners);
    corners     = [rows(1:n_keep), cols(1:n_keep), responses(1:n_keep)];
    corners_sub = corners(:, 1:2);
end
