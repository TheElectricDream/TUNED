function debugFeatureViewer(frame_curr, corners_curr, frameIndex, varargin)
%debugFeatureViewer Algorithm-agnostic interactive feature debug viewer.
%
%   debugFeatureViewer(FRAME_CURR, CORNERS_CURR, FRAMEINDEX)
%   debugFeatureViewer(..., 'Name', Value)
%
%   Provides two display modes for debugging any feature detection or
%   tracking algorithm:
%
%     'detect'  — Single frame with all detected features marked.
%                 Suitable for detection-only algorithms (e.g. Arc*).
%
%     'track'   — Side-by-side previous/current frames with coloured
%                 lines connecting matched features across frames.
%                 Suitable for detect+track pipelines (e.g. Harris-CV).
%
%   Execution blocks until the user presses "Next" or closes the
%   figure, unless 'auto_advance' is true.
%
%   Inputs (required):
%     frame_curr   - [nRows x nCols] current-frame image in SENSOR
%                    orientation (not transposed). uint8 or double.
%     corners_curr - [M x 2] current-frame corners [row, col] in
%                    sensor coordinates.
%     frameIndex   - Scalar frame number (for the title).
%
%   Name-Value Parameters (common):
%     'mode'           - 'detect' (default) or 'track'.
%     'algorithm_name' - Display name for the title bar and panel
%                        labels (default: 'Feature').
%     'gap'            - Pixel width of separator between panels
%                        (track mode only). Default: 4
%     'marker_size'    - Size of corner markers. Default: 8
%     'line_width'     - Width of match lines (track mode). Default: 1.5
%     'fig_tag'        - Figure tag for reuse across calls.
%                        Default: 'FeatureDebugViewer'
%     'colormap'       - Score colourmap: 'rg' (default), 'jet',
%                        or 'parula'.
%     'show_unmatched' - Show unmatched corners in grey.
%                        Default: true
%     'auto_advance'   - If true, do not block. Default: false
%
%   Name-Value Parameters (track mode only — ignored in detect mode):
%     'frame_prev'   - [nRows x nCols] previous-frame image. If empty
%                      or all-zero, a blank frame is shown.
%     'corners_prev' - [N x 2] previous-frame corners [row, col].
%     'matches'      - [P x 2] index pairs:
%                      matches(k,1) -> corners_curr index,
%                      matches(k,2) -> corners_prev index.
%     'scores'       - [P x 1] match quality scores in [0, 1].
%     'score_label'  - Name of the score metric for the legend
%                      (default: 'Score').
%
%   Usage (detect mode — single frame, detection only):
%
%     plot.debugFeatureViewer( ...
%         normalized_output_frame, ...
%         corner_positions, ...
%         frameIndex, ...
%         'mode',           'detect', ...
%         'algorithm_name', 'Arc*');
%
%   Usage (track mode — two frames, detection + matching):
%
%     plot.debugFeatureViewer( ...
%         normalized_output_frame, ...
%         corners_curr_valid, ...
%         frameIndex, ...
%         'mode',           'track', ...
%         'algorithm_name', 'Harris-CV', ...
%         'frame_prev',     debug_prev_frame, ...
%         'corners_prev',   corners_prev_frame_for_display, ...
%         'matches',        match_pairs, ...
%         'scores',         match_scores, ...
%         'score_label',    'NCC');
%
%   Controls:
%     [Next]  - Advance to the next frame (resumes execution).
%     [Stop]  - Stop debugging for all remaining frames.
%     Close   - Closing the figure also resumes execution. A new
%               figure is created on the next call.
%
%   See also: features.detectArcStarCorners, features.detectHarrisCV,
%             features.matchCVDescriptors

    % ================================================================
    % 0. Parse inputs
    % ================================================================
    ip = inputParser;
    addRequired(ip, 'frame_curr');
    addRequired(ip, 'corners_curr');
    addRequired(ip, 'frameIndex');

    % Common parameters
    addParameter(ip, 'mode',           'detect',              @ischar);
    addParameter(ip, 'algorithm_name', 'Feature',             @ischar);
    addParameter(ip, 'gap',            4,                     @isscalar);
    addParameter(ip, 'marker_size',    8,                     @isscalar);
    addParameter(ip, 'line_width',     1.5,                   @isscalar);
    addParameter(ip, 'fig_tag',        'FeatureDebugViewer',  @ischar);
    addParameter(ip, 'colormap',       'rg',                  @ischar);
    addParameter(ip, 'show_unmatched', true,                  @islogical);
    addParameter(ip, 'auto_advance',   false,                 @islogical);

    % Track-mode parameters (ignored in detect mode)
    addParameter(ip, 'frame_prev',     [],                    @isnumeric);
    addParameter(ip, 'corners_prev',   [],                    @isnumeric);
    addParameter(ip, 'matches',        [],                    @isnumeric);
    addParameter(ip, 'scores',         [],                    @isnumeric);
    addParameter(ip, 'score_label',    'Score',               @ischar);

    parse(ip, frame_curr, corners_curr, frameIndex, varargin{:});

    mode           = lower(ip.Results.mode);
    alg_name       = ip.Results.algorithm_name;
    gap_width      = ip.Results.gap;
    marker_sz      = ip.Results.marker_size;
    line_w         = ip.Results.line_width;
    fig_tag        = ip.Results.fig_tag;
    cmap_choice    = ip.Results.colormap;
    show_unmatched = ip.Results.show_unmatched;
    auto_advance   = ip.Results.auto_advance;
    frame_prev     = ip.Results.frame_prev;
    corners_prev   = ip.Results.corners_prev;
    matches        = ip.Results.matches;
    scores         = ip.Results.scores;
    score_label    = ip.Results.score_label;

    % Validate mode string
    assert(ismember(mode, {'detect', 'track'}), ...
        'debugFeatureViewer:badMode', ...
        'mode must be ''detect'' or ''track''.');

    % ================================================================
    % 1. Persistent stop flag
    % ================================================================
    persistent stop_all_flag;
    if isempty(stop_all_flag)
        stop_all_flag = false;
    end

    if stop_all_flag
        return;
    end

    % ================================================================
    % 2. Prepare frame data
    % ================================================================
    % Frames arrive in SENSOR orientation [imgSz(1) x imgSz(2)] =
    % [nRows_sensor x nCols_sensor] (e.g. [640 x 480], portrait).
    % Transpose to DISPLAY orientation [nCols_sensor x nRows_sensor]
    % (e.g. [480 x 640], landscape) so imagesc renders correctly.
    % Corner coordinates [sensor_row, sensor_col] are swapped to
    % [display_row, display_col] = [sensor_col, sensor_row] so that
    % the existing plot(corner(:,2), corner(:,1)) mapping still
    % produces correct (x, y) = (display_col, display_row) output.
    frame_curr = double(frame_curr)';
    [nRows, nCols] = size(frame_curr);

    if ~isempty(corners_curr)
        corners_curr = corners_curr(:, [2 1]);
    end
    if ~isempty(corners_prev)
        corners_prev = corners_prev(:, [2 1]);
    end

    if strcmp(mode, 'track')
        % Two-panel layout: [prev | gap | curr]
        if isempty(frame_prev) || all(frame_prev(:) == 0)
            frame_prev = zeros(nRows, nCols);
        else
            frame_prev = double(frame_prev)';
        end

        mx = max(max(frame_prev(:)), max(frame_curr(:)));
        if mx > 0
            frame_prev_n = frame_prev / mx;
            frame_curr_n = frame_curr / mx;
        else
            frame_prev_n = frame_prev;
            frame_curr_n = frame_curr;
        end

        separator = ones(nRows, gap_width) * 0.3;
        composite = [frame_prev_n, separator, frame_curr_n];
        x_offset  = nCols + gap_width;
    else
        % Single-panel layout
        mx = max(frame_curr(:));
        if mx > 0
            frame_curr_n = frame_curr / mx;
        else
            frame_curr_n = frame_curr;
        end
        composite = frame_curr_n;
        x_offset  = 0;
    end

    % ================================================================
    % 3. Create or reuse figure
    % ================================================================
    hFig = findobj('Type', 'figure', 'Tag', fig_tag);
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', [alg_name ' Debug Viewer'], ...
            'Tag', fig_tag, ...
            'NumberTitle', 'off', ...
            'MenuBar', 'none', ...
            'ToolBar', 'figure', ...
            'Units', 'normalized', ...
            'Position', [0.05 0.15 0.9 0.7], ...
            'Color', [0.15 0.15 0.15]);
    else
        set(hFig, 'Name', [alg_name ' Debug Viewer']);
        figure(hFig);
    end

    clf(hFig);

    % ================================================================
    % 4. Display composite image
    % ================================================================
    hAx = axes('Parent', hFig, ...
        'Position', [0.02 0.10 0.96 0.82]);

    imagesc(hAx, composite);
    colormap(hAx, gray(256));
    axis(hAx, 'image');
    hold(hAx, 'on');

    set(hAx, 'XTick', [], 'YTick', [], 'Box', 'on', ...
        'XColor', [0.4 0.4 0.4], 'YColor', [0.4 0.4 0.4]);

    % ================================================================
    % 5. Mode-specific rendering
    % ================================================================
    if strcmp(mode, 'detect')
        renderDetectMode(hAx, corners_curr, nRows, nCols, ...
            marker_sz, alg_name, frameIndex);
    else
        renderTrackMode(hAx, corners_curr, corners_prev, ...
            matches, scores, nRows, nCols, x_offset, ...
            marker_sz, line_w, cmap_choice, show_unmatched, ...
            alg_name, score_label, frameIndex, composite);
    end

    hold(hAx, 'off');

    % ================================================================
    % 6. UI controls
    % ================================================================
    if ~auto_advance
        btn_data.advance = false;
        btn_data.stop    = false;
        guidata(hFig, btn_data);

        % "Next" button
        uicontrol('Parent', hFig, ...
            'Style', 'pushbutton', ...
            'String', 'Next  ▶', ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'Units', 'normalized', ...
            'Position', [0.42 0.005 0.12 0.055], ...
            'BackgroundColor', [0.2 0.5 0.2], ...
            'ForegroundColor', 'w', ...
            'Callback', @(~,~) onNext(hFig));

        % "Stop" button
        uicontrol('Parent', hFig, ...
            'Style', 'pushbutton', ...
            'String', 'Stop All', ...
            'FontSize', 10, ...
            'Units', 'normalized', ...
            'Position', [0.56 0.005 0.10 0.055], ...
            'BackgroundColor', [0.5 0.15 0.15], ...
            'ForegroundColor', 'w', ...
            'Callback', @(~,~) onStop(hFig));

        % Frame counter — adapts label to mode
        if strcmp(mode, 'track')
            frame_label = sprintf('Frame %d → %d', frameIndex-1, frameIndex);
        else
            frame_label = sprintf('Frame %d', frameIndex);
        end

        uicontrol('Parent', hFig, ...
            'Style', 'text', ...
            'String', frame_label, ...
            'FontSize', 11, ...
            'FontWeight', 'bold', ...
            'Units', 'normalized', ...
            'Position', [0.30 0.005 0.11 0.055], ...
            'BackgroundColor', [0.15 0.15 0.15], ...
            'ForegroundColor', [0.8 0.8 0.8], ...
            'HorizontalAlignment', 'center');

        drawnow;

        % Block until button press or figure close
        try
            uiwait(hFig);
        catch
            % Figure was deleted externally — proceed
        end

        % Check which button was pressed
        if isvalid(hFig)
            btn_data = guidata(hFig);
            if btn_data.stop
                stop_all_flag = true;
            end
        end
    else
        drawnow;
    end
end


% =====================================================================
%  DETECT MODE RENDERER
% =====================================================================
function renderDetectMode(hAx, corners_curr, nRows, nCols, ...
    marker_sz, alg_name, frameIndex)
%renderDetectMode Draw all detected features on a single frame.

    % Panel title
    text(hAx, nCols/2, -8, ...
        sprintf('%s — Frame %d  (%d features)', ...
        alg_name, frameIndex, size(corners_curr, 1)), ...
        'Color', [0.3 1.0 0.5], 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % Plot all corners
    if ~isempty(corners_curr)
        plot(hAx, ...
            corners_curr(:, 2), ...     % x = col
            corners_curr(:, 1), ...     % y = row
            '+', 'Color', [0 0.9 0.4], ...
            'MarkerSize', marker_sz, ...
            'LineWidth', 1.2);
    else
        text(hAx, nCols/2, nRows/2, ...
            'No features detected', ...
            'Color', [1 0.4 0.4], 'FontSize', 14, ...
            'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold');
    end

    % Feature count annotation
    text(hAx, 5, nRows + 12, ...
        sprintf('%d features detected', size(corners_curr, 1)), ...
        'Color', [0.3 1.0 0.5], 'FontSize', 9);
end


% =====================================================================
%  TRACK MODE RENDERER
% =====================================================================
function renderTrackMode(hAx, corners_curr, corners_prev, ...
    matches, scores, nRows, nCols, x_offset, ...
    marker_sz, line_w, cmap_choice, show_unmatched, ...
    alg_name, score_label, frameIndex, composite)
%renderTrackMode Draw matched features across two side-by-side frames.

    % Panel labels
    text(hAx, nCols/2, -8, ...
        sprintf('%s — Frame %d (previous)', alg_name, frameIndex-1), ...
        'Color', [0.7 0.7 1.0], 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(hAx, x_offset + nCols/2, -8, ...
        sprintf('%s — Frame %d (current)', alg_name, frameIndex), ...
        'Color', [0.3 1.0 0.5], 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % Separator line
    xline(hAx, nCols + (x_offset - nCols)/2, ...
        'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1, 'LineStyle', '--');

    % --- Identify matched vs unmatched corners ---
    n_matches = size(matches, 1);

    if ~isempty(corners_curr) && n_matches > 0
        matched_curr_idx = matches(:, 1);
    else
        matched_curr_idx = [];
    end
    if ~isempty(corners_prev) && n_matches > 0
        matched_prev_idx = matches(:, 2);
    else
        matched_prev_idx = [];
    end

    % --- Unmatched corners (grey) ---
    if show_unmatched
        if ~isempty(corners_prev)
            unmatched_prev = setdiff(1:size(corners_prev,1), matched_prev_idx);
            if ~isempty(unmatched_prev)
                plot(hAx, ...
                    corners_prev(unmatched_prev, 2), ...
                    corners_prev(unmatched_prev, 1), ...
                    '+', 'Color', [0.45 0.45 0.45], ...
                    'MarkerSize', marker_sz * 0.7, ...
                    'LineWidth', 0.8);
            end
        end

        if ~isempty(corners_curr)
            unmatched_curr = setdiff(1:size(corners_curr,1), matched_curr_idx);
            if ~isempty(unmatched_curr)
                plot(hAx, ...
                    corners_curr(unmatched_curr, 2) + x_offset, ...
                    corners_curr(unmatched_curr, 1), ...
                    '+', 'Color', [0.45 0.45 0.45], ...
                    'MarkerSize', marker_sz * 0.7, ...
                    'LineWidth', 0.8);
            end
        end
    end

    % --- Match lines with score-based colouring ---
    if n_matches > 0
        match_colors = scoreToColor(scores, cmap_choice);

        for k = 1:n_matches
            ci = matches(k, 1);
            pi = matches(k, 2);

            x_prev = corners_prev(pi, 2);
            y_prev = corners_prev(pi, 1);
            x_curr = corners_curr(ci, 2) + x_offset;
            y_curr = corners_curr(ci, 1);

            clr = match_colors(k, :);

            line(hAx, [x_prev, x_curr], [y_prev, y_curr], ...
                'Color', [clr, 0.7], ...
                'LineWidth', line_w, ...
                'LineStyle', '-');

            plot(hAx, x_prev, y_prev, 'o', ...
                'MarkerEdgeColor', clr, ...
                'MarkerSize', marker_sz, ...
                'LineWidth', 1.5);
            plot(hAx, x_curr, y_curr, 's', ...
                'MarkerEdgeColor', clr, ...
                'MarkerSize', marker_sz, ...
                'LineWidth', 1.5);
        end

        % Score legend
        if ~isempty(scores)
            text(hAx, 5, nRows - 5, ...
                sprintf('%d matches  |  %s: %.3f \\pm %.3f  [%.3f, %.3f]', ...
                n_matches, score_label, mean(scores), std(scores), ...
                min(scores), max(scores)), ...
                'Color', 'w', 'FontSize', 10, ...
                'VerticalAlignment', 'bottom', ...
                'BackgroundColor', [0 0 0 0.6], ...
                'Margin', 3);
        end
    else
        text(hAx, size(composite,2)/2, nRows/2, ...
            'No matches this frame', ...
            'Color', [1 0.4 0.4], 'FontSize', 14, ...
            'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold');
    end

    % Corner count annotations
    n_prev = size(corners_prev, 1);
    n_curr = size(corners_curr, 1);
    text(hAx, 5, nRows + 12, ...
        sprintf('Prev: %d features', n_prev), ...
        'Color', [0.7 0.7 1.0], 'FontSize', 9);
    text(hAx, x_offset + 5, nRows + 12, ...
        sprintf('Curr: %d features', n_curr), ...
        'Color', [0.3 1.0 0.5], 'FontSize', 9);
end


% =====================================================================
%  Callback: "Next" button
% =====================================================================
function onNext(hFig)
    if isvalid(hFig)
        btn_data = guidata(hFig);
        btn_data.advance = true;
        guidata(hFig, btn_data);
        uiresume(hFig);
    end
end


% =====================================================================
%  Callback: "Stop All" button
% =====================================================================
function onStop(hFig)
    if isvalid(hFig)
        btn_data = guidata(hFig);
        btn_data.stop = true;
        guidata(hFig, btn_data);
        uiresume(hFig);
    end
end


% =====================================================================
%  Helper: map scores in [0, 1] to RGB colours
% =====================================================================
function colors = scoreToColor(scores, cmap_choice)
%scoreToColor Map quality scores in [0, 1] to RGB colours.
%
%   'rg'     - Red (low) → Yellow (mid) → Green (high)
%   'jet'    - Jet colourmap
%   'parula' - Parula colourmap

    n = numel(scores);
    colors = zeros(n, 3);

    % Clamp to [0, 1]
    s = max(0, min(1, scores));

    switch lower(cmap_choice)
        case 'rg'
            for i = 1:n
                if s(i) < 0.5
                    t = s(i) * 2;
                    colors(i,:) = [1.0, t, 0.0];
                else
                    t = (s(i) - 0.5) * 2;
                    colors(i,:) = [1.0 - t, 1.0, 0.0];
                end
            end

        case 'jet'
            cmap = jet(256);
            idx = max(1, min(256, round(s * 255) + 1));
            colors = cmap(idx, :);

        case 'parula'
            cmap = parula(256);
            idx = max(1, min(256, round(s * 255) + 1));
            colors = cmap(idx, :);

        otherwise
            colors = repmat([0 0.8 0.8], n, 1);
    end
end