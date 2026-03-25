function [auto_threshold, diagnostics] = findElbowThreshold(score_map, varargin)
% FINDELBOWTHRESHOLD  Data-driven threshold via geometric elbow detection.
%
%   AUTO_THRESHOLD = FINDELBOWTHRESHOLD(SCORE_MAP) sweeps candidate
%   thresholds over SCORE_MAP, computes the mean local binary
%   variance of the surviving pixel mask at each threshold, and
%   returns the threshold at the geometric elbow of the resulting
%   curve. The elbow is defined as the point of maximum perpendicular
%   distance from the chord connecting the curve endpoints after
%   normalizing both axes to [0, 1], following the Kneedle
%   convention [1].
%
%   The local binary variance at each pixel is p(1-p), where p is
%   the fraction of active pixels within a local window. This is
%   computed via a single conv2 call per threshold — no KD-trees or
%   nearest-neighbour searches. It captures local spatial clustering:
%   mixed noise-and-signal regions produce high p(1-p) while
%   homogeneous regions (pure cluster or pure background) produce
%   values near zero. The mean across the image summarizes the
%   overall degree of spatial mixing at each threshold level.
%
%   Dense-scene handling: LBV = p(1-p) is symmetric about p = 0.5.
%   When the initial fill fraction is high (dense scene), the LBV
%   curve rises from near-zero (p ≈ 1) to a peak at ~50% pixel
%   survival, then falls. Only the descending portion carries the
%   noise-to-signal transition information. The function detects
%   the LBV peak and truncates the curve to the descending portion
%   before applying Kneedle, ensuring monotonicity regardless of
%   initial scene density.
%
%   [AUTO_THRESHOLD, DIAGNOSTICS] = FINDELBOWTHRESHOLD(SCORE_MAP)
%   also returns a struct for inspection and plotting.
%
%   [...] = FINDELBOWTHRESHOLD(SCORE_MAP, 'Name', Value) accepts
%   optional name-value arguments:
%     'NumThresholds' - (50) Number of candidate thresholds.
%     'WinHalfSize'   - (5) Local window half-size in pixels,
%                        giving an (2*WinHalfSize+1)^2 window.
%     'UpperBound'    - (false) If true, also compute an upper
%                        threshold to reject hot pixels (values
%                        that are too persistent). The upper pass
%                        operates only on pixels surviving the
%                        lower threshold.
%
%   Inputs:
%     score_map - [H x W] 2D map of per-pixel scores (single or
%                 double). Zero and NaN entries are inactive.
%
%   Outputs:
%     auto_threshold - Scalar lower threshold at the geometric
%                      elbow. If 'UpperBound' is true, this is a
%                      [1 x 2] vector: [th_lo, th_hi].
%     diagnostics    - Struct with fields:
%       .th_vec       - [1 x M] Candidate thresholds (ascending,
%                       after truncation).
%       .mean_lbv     - [1 x M] Mean local binary variance (after
%                       truncation).
%       .th_norm      - [1 x M] Thresholds normalized to [0, 1].
%       .lbv_norm     - [1 x M] LBV normalized to [0, 1].
%       .perp_dist    - [1 x M] Perpendicular distance from chord.
%       .elbow_idx    - Index into the vectors above.
%       .peak_idx     - Index of the LBV peak in the *full* sweep
%                       (before truncation). peak_idx == 1 means
%                       the curve was already monotonically
%                       decreasing (sparse scene).
%       .full_th_vec  - [1 x N_TH] Full threshold sweep (before
%                       truncation), for diagnostic plotting.
%       .full_mean_lbv- [1 x N_TH] Full LBV curve (before
%                       truncation), for diagnostic plotting.
%       .upper        - (Only if 'UpperBound' is true) Struct with
%                       the same fields, computed on the flipped
%                       score axis for the upper pass.
%
%   Algorithm:
%     1. Extract unique nonzero score values and build a linearly
%        spaced sweep vector from min to max.
%     2. At each candidate threshold, binarize the map (>= th).
%     3. Compute local density p via box-filter convolution:
%          p = conv2(mask, kernel, 'same') / numel(kernel)
%     4. Mean local binary variance = mean(p .* (1 - p)) over
%        the full image.
%     5. Find the LBV peak. If peak_idx > 1 (dense scene), truncate
%        the curve to the descending portion [peak_idx : end].
%        If peak_idx == 1, the curve is already monotonic.
%     6. Prune flat trailing regions (LBV == 0).
%     7. Normalize both axes to [0, 1] and find the elbow
%        (max perpendicular distance from chord).
%
%   Computational cost:
%     O(N_TH x H x W) via conv2. For 640 x 480 with N_TH = 50:
%     ~50 convolutions, typically < 10 ms total in MATLAB.
%     Upper bound adds another ~50 convolutions on sparser data.
%
%   References:
%     [1] V. Satopaa, J. Albrecht, D. Irwin, and B. Raghavan,
%         "Finding a 'Kneedle' in a Haystack: Detecting Knee Points
%         in System Behavior," Proc. 31st Int. Conf. Distributed
%         Computing Systems Workshops, pp. 166-171, 2011.
%         DOI: 10.1109/ICDCSW.2011.20
%
%   See also: conv2, coherence.computeCoherenceMask,
%             plot.plotElbowDiagnostics

    % ----------------------------------------------------------------
    % 0. Parse arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addRequired(p, 'score_map', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'NumThresholds', 50,    @(x) isscalar(x) && x >= 3);
    addParameter(p, 'WinHalfSize',   5,     @(x) isscalar(x) && x >= 1);
    addParameter(p, 'UpperBound',    false, @(x) islogical(x) || x == 0 || x == 1);
    parse(p, score_map, varargin{:});

    N_th      = p.Results.NumThresholds;
    win_hz    = p.Results.WinHalfSize;
    do_upper  = logical(p.Results.UpperBound);

    % ----------------------------------------------------------------
    % 1. Lower-bound elbow
    % ----------------------------------------------------------------
    [th_lo, diag_lo] = run_elbow(score_map, N_th, win_hz);

    % ----------------------------------------------------------------
    % 2. (Optional) Upper-bound elbow on surviving pixels
    % ----------------------------------------------------------------
    if do_upper && ~isnan(th_lo)
        upper_map = score_map;
        upper_map(score_map < th_lo) = 0;

        surviving = upper_map(upper_map > 0 & ~isnan(upper_map));

        if numel(surviving) >= 10
            max_val    = max(surviving);
            flipped    = upper_map;
            active     = upper_map > 0 & ~isnan(upper_map);
            flipped(active)  = max_val - upper_map(active);
            flipped(~active) = 0;

            [th_flip, diag_up] = run_elbow(flipped, N_th, win_hz);

            if ~isnan(th_flip)
                th_hi = max_val - th_flip;

                if th_hi <= th_lo
                    th_hi = NaN;
                    diag_up.degenerate = true;
                else
                    diag_up.degenerate = false;
                end

                diag_up.th_vec_original = max_val - fliplr(diag_up.th_vec);
                diag_up.mean_lbv_original = fliplr(diag_up.mean_lbv);
            else
                th_hi   = NaN;
                diag_up = struct('degenerate', true);
            end
        else
            th_hi   = NaN;
            diag_up = struct('degenerate', true);
        end

        auto_threshold = [th_lo, th_hi];

        if nargout > 1
            diagnostics       = diag_lo;
            diagnostics.upper = diag_up;
        end
    else
        auto_threshold = th_lo;
        if nargout > 1
            diagnostics = diag_lo;
        end
    end
end


% ====================================================================
%  HELPER: core elbow detection with dense-scene peak truncation
% ====================================================================
function [threshold, diag_out] = run_elbow(score_map, N_th, win_hz)

    vals = score_map(score_map > 0 & ~isnan(score_map));

    if numel(vals) < 10
        warning('findElbowThreshold:tooFewActive', ...
            'Score map has fewer than 10 nonzero pixels. Returning NaN.');
        threshold = NaN;
        diag_out  = struct();
        return;
    end

    th_vec = linspace(min(vals), max(vals), N_th);

    % Box-filter kernel
    win_side = 2 * win_hz + 1;
    kernel   = ones(win_side) / win_side^2;

    % Sweep: mean local binary variance at each threshold
    mean_lbv = zeros(1, N_th);
    for k = 1:N_th
        mask = double(score_map >= th_vec(k));
        px   = conv2(mask, kernel, 'same');
        lbv  = px .* (1 - px);
        mean_lbv(k) = mean(lbv(:));
    end

    % Store the full curve before any truncation
    full_th_vec   = th_vec;
    full_mean_lbv = mean_lbv;

    % ----------------------------------------------------------------
    % Dense-scene handling: truncate to descending portion.
    %
    % LBV = p(1-p) is symmetric about p = 0.5. In sparse scenes the
    % curve starts high and decreases monotonically (p moves from
    % moderate toward 0), so peak_idx == 1 and nothing is removed.
    % In dense scenes the curve rises first (p moves from ~1 toward
    % 0.5) then falls. Only the descending portion — from the peak
    % onward — captures the transition from spatially mixed
    % (noise + signal) to homogeneous clusters. We discard the
    % ascending portion so that Kneedle sees a monotonic curve.
    % ----------------------------------------------------------------
    [~, peak_idx] = max(mean_lbv);
    th_vec   = th_vec(peak_idx:end);
    mean_lbv = mean_lbv(peak_idx:end);

    % Prune flat trailing regions (LBV == 0)
    good = mean_lbv > 0;
    if sum(good) < 3
        threshold = th_vec(round(numel(th_vec) / 2));
        diag_out  = struct();
        diag_out.peak_idx      = peak_idx;
        diag_out.full_th_vec   = full_th_vec;
        diag_out.full_mean_lbv = full_mean_lbv;
        return;
    end

    th_g  = th_vec(good);
    lbv_g = mean_lbv(good);

    % Normalize both axes to [0, 1]
    th_range  = th_g(end) - th_g(1);
    lbv_range = max(lbv_g) - min(lbv_g);

    if th_range == 0 || lbv_range == 0
        threshold = th_g(round(numel(th_g) / 2));
        diag_out  = struct();
        diag_out.peak_idx      = peak_idx;
        diag_out.full_th_vec   = full_th_vec;
        diag_out.full_mean_lbv = full_mean_lbv;
        return;
    end

    th_n  = (th_g  - th_g(1))    / th_range;
    lbv_n = (lbv_g - min(lbv_g)) / lbv_range;

    % Geometric elbow: max perpendicular distance from chord
    v_chord = [th_n(end) - th_n(1), lbv_n(end) - lbv_n(1)];
    v_hat   = v_chord / norm(v_chord);

    w_x = th_n - th_n(1);
    w_y = lbv_n - lbv_n(1);
    perp_dist = abs(w_x * v_hat(2) - w_y * v_hat(1));

    [~, elbow_idx] = max(perp_dist);
    threshold = th_g(elbow_idx);

    % Pack diagnostics
    diag_out.th_vec        = th_g;
    diag_out.mean_lbv      = lbv_g;
    diag_out.th_norm       = th_n;
    diag_out.lbv_norm      = lbv_n;
    diag_out.perp_dist     = perp_dist;
    diag_out.elbow_idx     = elbow_idx;
    diag_out.peak_idx      = peak_idx;
    diag_out.full_th_vec   = full_th_vec;
    diag_out.full_mean_lbv = full_mean_lbv;
end