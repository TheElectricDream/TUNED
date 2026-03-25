function fig = plotElbowDiagnostics(auto_threshold, diagnostics, varargin)
% PLOTELBOWDIAGNOSTICS  Diagnostic figure for elbow threshold(s).
%
%   FIG = PLOTELBOWDIAGNOSTICS(AUTO_THRESHOLD, DIAGNOSTICS) produces a
%   diagnostic figure. The left panel shows the full LBV sweep with
%   the LBV peak (truncation point) and selected elbow marked. The
%   right panel shows the normalized descending curve used by Kneedle,
%   with the chord and perpendicular distance peak.
%
%   If AUTO_THRESHOLD is a [1 x 2] vector (dual-threshold mode) and
%   DIAGNOSTICS contains an .upper field, a second row of panels is
%   added for the upper-bound elbow detection.
%
%   FIG = PLOTELBOWDIAGNOSTICS(..., 'Name', Value) accepts optional
%   name-value arguments:
%     'Title'       - ('') String prepended to subplot titles.
%     'FigPosition' - ([]) Figure position vector. Default adapts
%                     to the number of panels.
%
%   Inputs:
%     auto_threshold - Scalar elbow threshold, or [1 x 2] vector
%                      [th_lo, th_hi] from dual-threshold mode.
%     diagnostics    - Struct returned by findElbowThreshold.
%
%   Outputs:
%     fig - Figure handle.
%
%   See also: stats.findElbowThreshold

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'Title', '', @ischar);
    addParameter(p, 'FigPosition', [], ...
        @(x) isempty(x) || (isnumeric(x) && numel(x) == 4));
    parse(p, varargin{:});
    label_prefix = p.Results.Title;
    fig_pos      = p.Results.FigPosition;

    has_upper = numel(auto_threshold) == 2 && ...
                isfield(diagnostics, 'upper') && ...
                isstruct(diagnostics.upper) && ...
                (~isfield(diagnostics.upper, 'degenerate') || ...
                 ~diagnostics.upper.degenerate);

    if ~isempty(label_prefix)
        label_prefix = [label_prefix ' — '];
    end

    % Determine layout
    if has_upper
        n_rows = 2;
        if isempty(fig_pos)
            fig_pos = [100 100 900 650];
        end
    else
        n_rows = 1;
        if isempty(fig_pos)
            fig_pos = [100 100 900 350];
        end
    end

    fig = figure('Position', fig_pos);

    % === Row 1: Lower-bound elbow ====================================
    th_lo = auto_threshold(1);

    % Unpack diagnostics
    th_g  = diagnostics.th_vec;
    lbv_g = diagnostics.mean_lbv;
    th_n  = diagnostics.th_norm;
    lbv_n = diagnostics.lbv_norm;
    ei    = diagnostics.elbow_idx;

    has_full = isfield(diagnostics, 'full_th_vec') && ...
               isfield(diagnostics, 'full_mean_lbv');
    has_peak = isfield(diagnostics, 'peak_idx');

    % --- Left panel: full curve with truncation + elbow marked ---
    subplot(n_rows, 2, 1);
    hold on;

    if has_full
        % Show the full sweep in light gray
        plot(diagnostics.full_th_vec, diagnostics.full_mean_lbv, ...
            '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, ...
            'DisplayName', 'Full sweep');

        % Mark the truncation point (LBV peak)
        if has_peak && diagnostics.peak_idx > 1
            pk = diagnostics.peak_idx;
            plot(diagnostics.full_th_vec(pk), ...
                diagnostics.full_mean_lbv(pk), ...
                'ks', 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 0.5], ...
                'DisplayName', sprintf('LBV peak (idx %d)', pk));
        end
    end

    % Overlay the descending portion used by Kneedle
    plot(th_g, lbv_g, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2, ...
        'DisplayName', 'Descending portion');
    xline(th_lo, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(th_lo, lbv_g(ei), 'rp', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'r', ...
        'DisplayName', sprintf('Elbow @ %.4f', th_lo));

    xlabel('Threshold');
    ylabel('Mean local binary variance');
    if has_peak && diagnostics.peak_idx > 1
        title(sprintf('%sLower — elbow @ %.4f (peak truncated)', ...
            label_prefix, th_lo));
    else
        title(sprintf('%sLower — elbow @ %.4f', label_prefix, th_lo));
    end
    legend('Location', 'best');
    grid on;
    hold off;

    % --- Right panel: normalized + chord ---
    subplot(n_rows, 2, 2);
    plot(th_n, lbv_n, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
    hold on;
    plot([th_n(1) th_n(end)], [lbv_n(1) lbv_n(end)], 'k--', ...
        'LineWidth', 1.0);
    plot(th_n(ei), lbv_n(ei), 'rp', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'r');
    draw_perp_line(th_n, lbv_n, ei);

    xlabel('Threshold (normalized)');
    ylabel('Mean LBV (normalized)');
    title(sprintf('%sLower — max \\perp dist', label_prefix));
    grid on;
    hold off;

    % === Row 2: Upper-bound elbow (if enabled) =======================
    if has_upper
        th_hi = auto_threshold(2);
        du    = diagnostics.upper;

        % Use original-axis remapped vectors for intuitive display
        th_g_u  = du.th_vec_original;
        lbv_g_u = du.mean_lbv_original;
        ei_u    = numel(du.th_vec) - du.elbow_idx + 1;

        subplot(n_rows, 2, 3);
        hold on;

        if isfield(du, 'full_th_vec') && isfield(du, 'full_mean_lbv')
            % Show full flipped sweep in gray (remapped to original axis)
            max_val = max(du.full_th_vec);
            plot(max_val - fliplr(du.full_th_vec), ...
                fliplr(du.full_mean_lbv), ...
                '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0);
        end

        plot(th_g_u, lbv_g_u, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
        xline(th_hi, 'm--', 'LineWidth', 1.5);
        xline(th_lo, 'r:', 'LineWidth', 1.0, 'Alpha', 0.4);
        plot(th_hi, lbv_g_u(ei_u), 'mp', ...
            'MarkerSize', 14, 'MarkerFaceColor', 'm');

        xlabel('Threshold (original axis)');
        ylabel('Mean local binary variance');
        title(sprintf('%sUpper — elbow @ %.4f', label_prefix, th_hi));
        legend('Full sweep', 'Descending', 'Upper', 'Lower', ...
            'Location', 'best');
        grid on;
        hold off;

        subplot(n_rows, 2, 4);
        plot(du.th_norm, du.lbv_norm, 'b-o', 'MarkerSize', 3, ...
            'LineWidth', 1.2);
        hold on;
        plot([du.th_norm(1) du.th_norm(end)], ...
             [du.lbv_norm(1) du.lbv_norm(end)], 'k--', 'LineWidth', 1.0);
        plot(du.th_norm(du.elbow_idx), du.lbv_norm(du.elbow_idx), ...
            'mp', 'MarkerSize', 14, 'MarkerFaceColor', 'm');
        draw_perp_line(du.th_norm, du.lbv_norm, du.elbow_idx);

        xlabel('Flipped threshold (normalized)');
        ylabel('Mean LBV (normalized)');
        title(sprintf('%sUpper — max \\perp dist (flipped)', label_prefix));
        grid on;
        hold off;
    end
end


% ====================================================================
%  HELPER: draw perpendicular line from elbow to chord
% ====================================================================
function draw_perp_line(th_n, lbv_n, ei)
    p1 = [th_n(1), lbv_n(1)];
    v  = [th_n(end), lbv_n(end)] - p1;
    w  = [th_n(ei), lbv_n(ei)] - p1;
    proj_len = dot(w, v) / dot(v, v);
    foot = p1 + proj_len * v;
    plot([th_n(ei) foot(1)], [lbv_n(ei) foot(2)], 'r-', ...
        'LineWidth', 1.0);
end