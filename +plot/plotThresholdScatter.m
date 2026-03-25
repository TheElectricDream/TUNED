function [fig, auto_threshold, diagnostics] = plotThresholdScatter(map, varargin)
% PLOTTHRESHOLDSCATTER  3D scatter plot with elbow-threshold plane.
%
%   FIG = PLOTTHRESHOLDSCATTER(MAP) converts a 2D value map to a
%   point cloud via process.generateMeshFromFrame, computes the
%   data-driven elbow threshold via stats.findElbowThreshold, and
%   displays a 3D scatter plot where points below and above the
%   threshold are rendered in distinct colors. A semi-transparent
%   horizontal plane at z = auto_threshold shows the decision
%   boundary.
%
%   [FIG, AUTO_THRESHOLD, DIAGNOSTICS] = PLOTTHRESHOLDSCATTER(MAP)
%   also returns the scalar threshold and the diagnostics struct
%   from stats.findElbowThreshold for downstream use or plotting
%   with plot.plotElbowDiagnostics.
%
%   [...] = PLOTTHRESHOLDSCATTER(MAP, 'Name', Value) accepts
%   optional name-value arguments:
%     'NThresholds'  - (50)  Number of candidate thresholds passed
%                      to stats.findElbowThreshold.
%     'WindowSize'   - (5)   Local window half-size [px] passed to
%                      stats.findElbowThreshold.
%     'MarkerSize'   - (36)  Scatter marker area [pt^2].
%     'PlaneAlpha'   - (0.25) Transparency of the threshold plane.
%     'PlaneColor'   - ([1 0 0]) RGB color of the threshold plane.
%     'BelowColor'   - ([0.2 0.4 0.8]) RGB for points below threshold.
%     'AboveColor'   - ([0.9 0.2 0.1]) RGB for points above threshold.
%     'Title'        - ('') String prepended to the figure title.
%     'FigPosition'  - ([100 100 800 600]) Figure position vector.
%
%   Inputs:
%     map - [H x W] 2D map of per-pixel scores (e.g., time surface,
%           IEI mean, CV map). Zero-valued pixels are treated as
%           inactive and excluded from the scatter.
%
%   Outputs:
%     fig            - Figure handle.
%     auto_threshold - Scalar elbow threshold.
%     diagnostics    - Struct from stats.findElbowThreshold.
%
%   Example:
%     % Visualize a time surface with its elbow threshold
%     [fig, th, diag] = plot.plotThresholdScatter(surface_map, ...
%         'Title', 'IEI-ATS Surface');
%     plot.plotElbowDiagnostics(th, diag, 'Title', 'IEI-ATS Surface');
%
%   See also: stats.findElbowThreshold, plot.mapToScatterPlot,
%             plot.plotElbowDiagnostics, process.generateMeshFromFrame

    % ----------------------------------------------------------------
    % 0. Parse arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addRequired(p, 'map', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'NThresholds', 50, ...
        @(x) isnumeric(x) && isscalar(x) && x > 2);
    addParameter(p, 'WindowSize', 5, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'MarkerSize', 36, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'PlaneAlpha', 0.25, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    addParameter(p, 'PlaneColor', [1 0 0], ...
        @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'BelowColor', [0.2 0.4 0.8], ...
        @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'AboveColor', [0.9 0.2 0.1], ...
        @(x) isnumeric(x) && numel(x) == 3);
    addParameter(p, 'Title', '', @ischar);
    addParameter(p, 'FigPosition', [100 100 800 600], ...
        @(x) isnumeric(x) && numel(x) == 4);
    parse(p, map, varargin{:});

    n_th       = p.Results.NThresholds;
    win_sz     = p.Results.WindowSize;
    mk_sz      = p.Results.MarkerSize;
    plane_a    = p.Results.PlaneAlpha;
    plane_c    = p.Results.PlaneColor;
    col_below  = p.Results.BelowColor;
    col_above  = p.Results.AboveColor;
    label_pre  = p.Results.Title;
    fig_pos    = p.Results.FigPosition;

    % ----------------------------------------------------------------
    % 1. Compute elbow threshold
    % ----------------------------------------------------------------
    [auto_threshold, diagnostics] = ...
        stats.findElbowThreshold(map, 'NumThresholds', n_th, 'WinHalfSize', win_sz);

    if isnan(auto_threshold)
        warning('plotThresholdScatter:nanThreshold', ...
            'Elbow threshold returned NaN — too few active pixels.');
        fig = figure('Position', fig_pos);
        text(0.5, 0.5, 'Insufficient active pixels for threshold.', ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
        axis off;
        return;
    end

    % ----------------------------------------------------------------
    % 2. Convert map to point cloud (reuse project utility)
    % ----------------------------------------------------------------
    pointCloud = process.generateMeshFromFrame(map);

    % Strip NaN rows (inactive pixels)
    valid = ~isnan(pointCloud(:, 3));
    x = pointCloud(valid, 1);
    y = pointCloud(valid, 2);
    z = pointCloud(valid, 3);

    % ----------------------------------------------------------------
    % 3. Partition into below / above threshold
    % ----------------------------------------------------------------
    below = z < auto_threshold;
    above = ~below;

    % ----------------------------------------------------------------
    % 4. Build figure
    % ----------------------------------------------------------------
    fig = figure('Position', fig_pos);
    hold on;

    % Points below threshold
    if any(below)
        scatter3(x(below), y(below), z(below), ...
            mk_sz, repmat(col_below, sum(below), 1), 'filled', ...
            'MarkerFaceAlpha', 0.6, ...
            'DisplayName', sprintf('Below (n=%d)', sum(below)));
    end

    % Points above threshold
    if any(above)
        scatter3(x(above), y(above), z(above), ...
            mk_sz, repmat(col_above, sum(above), 1), 'filled', ...
            'MarkerFaceAlpha', 0.8, ...
            'DisplayName', sprintf('Above (n=%d)', sum(above)));
    end

    % ----------------------------------------------------------------
    % 5. Threshold plane
    % ----------------------------------------------------------------
    x_lim = [min(x), max(x)];
    y_lim = [min(y), max(y)];

    patch_x = [x_lim(1) x_lim(2) x_lim(2) x_lim(1)];
    patch_y = [y_lim(1) y_lim(1) y_lim(2) y_lim(2)];
    patch_z = auto_threshold * ones(1, 4);

    patch(patch_x, patch_y, patch_z, plane_c, ...
        'FaceAlpha', plane_a, ...
        'EdgeColor', plane_c, ...
        'LineWidth', 1.0, ...
        'DisplayName', sprintf('Threshold = %.4f', auto_threshold));

    % ----------------------------------------------------------------
    % 6. Labels and formatting
    % ----------------------------------------------------------------
    if ~isempty(label_pre)
        label_pre = [label_pre ' — '];
    end

    xlabel('X [px]');
    ylabel('Y [px]');
    zlabel('Score');
    title(sprintf('%sElbow threshold @ %.4f', label_pre, auto_threshold));
    legend('Location', 'best');
    colormap(jet);
    view(3);
    grid on;
    set(gca, 'FontSize', 16);
    hold off;

end