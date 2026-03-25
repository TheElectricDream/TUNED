function fig = plotSelectedMetric(loopLog, filterNames, metricName, varargin)
% PLOTSELECTEDMETRIC  Time-series comparison of a single metric across filters.
%
%   FIG = PLOTSELECTEDMETRIC(LOOPLOG, FILTERNAMES, METRICNAME) extracts
%   the per-frame vector for METRICNAME from each filter's log and
%   overlays them on a single figure for visual comparison.
%
%   FIG = PLOTSELECTEDMETRIC(..., 'Name', Value) accepts optional
%   name-value arguments:
%     'Title'      - ('') String prepended to the plot title.
%     'LineWidth'  - (1.5) Line width for all traces.
%     'FigPosition' - ([100 100 800 450]) Figure position vector.
%
%   Inputs:
%     loopLog     - {K x 1} Cell array. Each cell contains a
%                   [1 x F] struct array of per-frame metrics
%                   (output of metrics.computeFrameMetrics).
%     filterNames - {K x 1} Cell array of strings. Display name
%                   for each filter (used in the legend).
%     metricName  - Char. Field name to extract from loopLog
%                   entries. Must match a field in the per-frame
%                   metrics struct 
%                   See METRICLABEL (local) for supported names.
%
%   Outputs:
%     fig - Figure handle.
%
%   Example:
%     fig = plot.plotSelectedMetric(loopLog, filterNames, 'srr');
%     fig = plot.plotSelectedMetric(loopLog, filterNames, ...
%               'contrast_iqr', 'Title', 'Nominal lighting');
%
%   See also: metrics.computeFrameMetrics, metrics.summarizeMetrics

    % ----------------------------------------------------------------
    % Parse optional arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addRequired(p, 'loopLog',     @iscell);
    addRequired(p, 'filterNames', @iscell);
    addRequired(p, 'metricName',  @ischar);
    addParameter(p, 'Title',       '', @ischar);
    addParameter(p, 'LineWidth',   1.5, @isnumeric);
    addParameter(p, 'FigPosition', [655         258        1355         851], ...
        @(x) isnumeric(x) && numel(x) == 4);
    parse(p, loopLog, filterNames, metricName, varargin{:});

    label_prefix = p.Results.Title;
    lw           = p.Results.LineWidth;
    fig_pos      = p.Results.FigPosition;

    % ----------------------------------------------------------------
    % Resolve display label and y-axis label for the metric
    % ----------------------------------------------------------------
    [y_label, display_name] = metricLabel(metricName);

    if ~isempty(label_prefix)
        plot_title = sprintf('%s — %s', label_prefix, display_name);
    else
        plot_title = display_name;
    end

    % ----------------------------------------------------------------
    % Plot
    % ----------------------------------------------------------------
    fig = figure('Position', fig_pos);
    hold on;

    for k = 1:numel(filterNames)
        vec_k = [loopLog{k}.(metricName)];
        plot(vec_k, 'DisplayName', filterNames{k}, 'LineWidth', lw);
        axis tight;
        grid on;
    end

    xlabel('Frame Index');
    ylabel(y_label);
    title(plot_title);
    legend('Location', 'best');
    grid on;
    set(gca, 'FontSize', 16);
    hold off;

end

% ====================================================================
% Local function: metric metadata lookup
% ====================================================================
function [y_label, display_name] = metricLabel(name)
% METRICLABEL  Return axis label and display name for a metric field.
%
%   Centralises all per-metric display strings so that adding a new
%   metric only requires one new entry here.

    lookup = { ...
        % field name                    y-axis label                        display name
        'SRR',                          'Signal Retention Rate',            'Signal Retention Rate (SRR)'; 
        'ClarkEvansRemaining',                   'Clark-Evans (Background)',                      'Clark Evans';
        'ClarkEvansRemoved',                   'Clark-Evans (Foreground)',                      'Clark Evans';
        'ComputeTimeFilter',            'Compute Time (s)',                 'Compute Time (Filter)';
        'ComputeTimeAccumulator',       'Compute Time (s)',                 'Compute Time (Accumulator)';
        'EventsInFrame',                'Event Count',                      'Total Events per Frame';
        'FilteredEvents',               'Filtered Event Count',             'Total Filtered Events per Frame';
        'FilteringMEVs',                'MEVs / s',                         'Filtering MEVs (events per second)';
        'AccumulatorMEVs',              'MEVs / s',                         'Accumulator MEVs (events per second)';
    };

    idx = find(strcmp(lookup(:,1), name), 1);
    if ~isempty(idx)
        y_label      = lookup{idx, 2};
        display_name = lookup{idx, 3};
    else
        warning('plotSelectedMetric:unknownMetric', ...
            'Unrecognized metric "%s". Using raw field name.', name);
        y_label      = name;
        display_name = name;
    end

end