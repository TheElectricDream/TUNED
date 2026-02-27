function debugMatchViewer(frame_prev, frame_curr, ...
    corners_prev, corners_curr, matches, scores, frameIndex, varargin)
%debugMatchViewer Interactive side-by-side match visualisation.
%
%   debugMatchViewer(FRAME_PREV, FRAME_CURR, CORNERS_PREV, CORNERS_CURR,
%       MATCHES, SCORES, FRAMEINDEX)
%   debugMatchViewer(..., 'Name', Value)
%
%   Displays two consecutive frames side-by-side with coloured lines
%   connecting matched features across frames. Execution blocks until
%   the user presses the "Next" button or closes the figure.
%
%   The left panel shows frame N-1 with previous-frame corners marked.
%   The right panel shows frame N with current-frame corners marked.
%   Lines span from matched corners in the left panel to their
%   correspondences in the right panel, coloured by NCC score
%   (red = weak, green = strong).
%
%   Inputs:
%     frame_prev   - [nRows x nCols] previous-frame image in SENSOR
%                    orientation (not transposed). Can be uint8 or
%                    double. If empty or all-zero, a blank frame is
%                    shown (useful for frameIndex == 1).
%     frame_curr   - [nRows x nCols] current-frame image in SENSOR
%                    orientation.
%     corners_prev - [N x 2] previous-frame corners [row, col] in
%                    sensor coordinates (from extractCVDescriptors).
%     corners_curr - [M x 2] current-frame corners [row, col] in
%                    sensor coordinates.
%     matches      - [P x 2] index pairs from matchCVDescriptors:
%                    matches(k,1) -> corners_curr index,
%                    matches(k,2) -> corners_prev index.
%     scores       - [P x 1] NCC scores for each match.
%     frameIndex   - Scalar frame number (for the title).
%
%   Name-Value Parameters:
%     'gap'          - Pixel width of the separator between panels.
%                      Default: 4
%     'marker_size'  - Size of corner markers. Default: 8
%     'line_width'   - Width of match lines. Default: 1.5
%     'fig_tag'      - Figure tag for reuse across calls. If a figure
%                      with this tag exists, it is reused (position
%                      preserved). Default: 'HarrisCVDebugViewer'
%     'colormap'     - Colormap for NCC score. 'rg' for red-green
%                      (default), 'jet', or 'parula'.
%     'show_unmatched' - Show unmatched corners with grey markers.
%                        Default: true
%     'auto_advance' - If true, do not block (useful for batch
%                      recording). Default: false
%
%   Usage in main.m (inside the main loop, after matching):
%
%     plot.debugMatchViewer( ...
%         frame_prev_for_debug, ...
%         normalized_output_frame, ...
%         corners_prev_frame_for_display, ...
%         corners_curr_valid, ...
%         harris_cv_matches, ...
%         harris_cv_scores, ...
%         frameIndex);
%
%   Controls:
%     [Next]  - Advance to the next frame (resumes execution).
%     [Stop]  - Stop debugging for all remaining frames.
%     Close   - Closing the figure also resumes execution. A new
%               figure is created on the next call.
%
%   See also: features.matchCVDescriptors, features.detectHarrisCV

    % ----------------------------------------------------------------
    % 0. Parse inputs
    % ----------------------------------------------------------------
    ip = inputParser;
    addRequired(ip, 'frame_prev');
    addRequired(ip, 'frame_curr');
    addRequired(ip, 'corners_prev');
    addRequired(ip, 'corners_curr');
    addRequired(ip, 'matches');
    addRequired(ip, 'scores');
    addRequired(ip, 'frameIndex');
    addParameter(ip, 'gap',            4,                      @isscalar);
    addParameter(ip, 'marker_size',    8,                      @isscalar);
    addParameter(ip, 'line_width',     1.5,                    @isscalar);
    addParameter(ip, 'fig_tag',        'HarrisCVDebugViewer',  @ischar);
    addParameter(ip, 'colormap',       'rg',                   @ischar);
    addParameter(ip, 'show_unmatched', true,                   @islogical);
    addParameter(ip, 'auto_advance',   false,                  @islogical);
    parse(ip, frame_prev, frame_curr, corners_prev, ...
        corners_curr, matches, scores, frameIndex, varargin{:});

    gap_width      = ip.Results.gap;
    marker_sz      = ip.Results.marker_size;
    line_w         = ip.Results.line_width;
    fig_tag        = ip.Results.fig_tag;
    cmap_choice    = ip.Results.colormap;
    show_unmatched = ip.Results.show_unmatched;
    auto_advance   = ip.Results.auto_advance;

    % ----------------------------------------------------------------
    % 1. Prepare frames
    % ----------------------------------------------------------------
    frame_curr = double(frame_curr);
    [nRows, nCols] = size(frame_curr);

    if isempty(frame_prev) || all(frame_prev(:) == 0)
        frame_prev = zeros(nRows, nCols);
    else
        frame_prev = double(frame_prev);
    end

    % Normalise both to [0, 1] for display
    mx = max(max(frame_prev(:)), max(frame_curr(:)));
    if mx > 0
        frame_prev_n = frame_prev / mx;
        frame_curr_n = frame_curr / mx;
    else
        frame_prev_n = frame_prev;
        frame_curr_n = frame_curr;
    end

    % Build composite: [prev | gap | curr]
    separator = ones(nRows, gap_width) * 0.3;   % dark grey bar
    composite = [frame_prev_n, separator, frame_curr_n];

    % x-offset for the right panel
    x_offset = nCols + gap_width;

    % ----------------------------------------------------------------
    % 2. Create or reuse figure
    % ----------------------------------------------------------------
    persistent stop_all_flag;
    if isempty(stop_all_flag)
        stop_all_flag = false;
    end

    % Check stop flag from previous "Stop" press
    if stop_all_flag
        return;
    end

    hFig = findobj('Type', 'figure', 'Tag', fig_tag);
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Harris-CV Match Debugger', ...
            'Tag', fig_tag, ...
            'NumberTitle', 'off', ...
            'MenuBar', 'none', ...
            'ToolBar', 'figure', ...
            'Units', 'normalized', ...
            'Position', [0.05 0.15 0.9 0.7], ...
            'Color', [0.15 0.15 0.15]);
    else
        figure(hFig);
    end

    clf(hFig);

    % ----------------------------------------------------------------
    % 3. Display composite image
    % ----------------------------------------------------------------
    hAx = axes('Parent', hFig, ...
        'Position', [0.02 0.10 0.96 0.82]);

    % imagesc uses (x, y) = (col, row), so for sensor-space images:
    %   x-axis = columns (1 to nCols for left, offset for right)
    %   y-axis = rows (1 to nRows, top-to-bottom)
    imagesc(hAx, composite);
    colormap(hAx, gray(256));
    axis(hAx, 'image');
    hold(hAx, 'on');

    % Panel labels
    text(hAx, nCols/2, -8, sprintf('Frame %d (previous)', frameIndex-1), ...
        'Color', [0.7 0.7 1.0], 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(hAx, x_offset + nCols/2, -8, sprintf('Frame %d (current)', frameIndex), ...
        'Color', [0.3 1.0 0.5], 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % Separator label
    xline(hAx, nCols + gap_width/2, 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1, 'LineStyle', '--');

    set(hAx, 'XTick', [], 'YTick', [], 'Box', 'on', ...
        'XColor', [0.4 0.4 0.4], 'YColor', [0.4 0.4 0.4]);

    % ----------------------------------------------------------------
    % 4. Identify matched vs unmatched corners
    % ----------------------------------------------------------------
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

    % ----------------------------------------------------------------
    % 5. Plot unmatched corners (grey)
    % ----------------------------------------------------------------
    if show_unmatched
        % Unmatched previous-frame corners
        if ~isempty(corners_prev)
            unmatched_prev = setdiff(1:size(corners_prev,1), matched_prev_idx);
            if ~isempty(unmatched_prev)
                plot(hAx, ...
                    corners_prev(unmatched_prev, 2), ...    % x = col
                    corners_prev(unmatched_prev, 1), ...    % y = row
                    '+', 'Color', [0.45 0.45 0.45], ...
                    'MarkerSize', marker_sz * 0.7, ...
                    'LineWidth', 0.8);
            end
        end

        % Unmatched current-frame corners (offset to right panel)
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

    % ----------------------------------------------------------------
    % 6. Draw match lines with NCC-based colouring
    % ----------------------------------------------------------------
    if n_matches > 0
        % Build colour array
        match_colors = nccToColor(scores, cmap_choice);

        for k = 1:n_matches
            ci = matches(k, 1);   % index into corners_curr
            pi = matches(k, 2);   % index into corners_prev

            % Previous corner position (left panel)
            x_prev = corners_prev(pi, 2);           % col
            y_prev = corners_prev(pi, 1);            % row

            % Current corner position (right panel, offset)
            x_curr = corners_curr(ci, 2) + x_offset;
            y_curr = corners_curr(ci, 1);

            clr = match_colors(k, :);

            % Draw connecting line
            line(hAx, [x_prev, x_curr], [y_prev, y_curr], ...
                'Color', [clr, 0.7], ...    % with alpha via 4th element
                'LineWidth', line_w, ...
                'LineStyle', '-');

            % Matched corner markers (coloured)
            plot(hAx, x_prev, y_prev, 'o', ...
                'MarkerEdgeColor', clr, ...
                'MarkerSize', marker_sz, ...
                'LineWidth', 1.5);
            plot(hAx, x_curr, y_curr, 's', ...
                'MarkerEdgeColor', clr, ...
                'MarkerSize', marker_sz, ...
                'LineWidth', 1.5);
        end

        % NCC score legend (text in bottom-left)
        if ~isempty(scores)
            text(hAx, 5, nRows - 5, ...
                sprintf('%d matches  |  NCC: %.3f \\pm %.3f  [%.3f, %.3f]', ...
                n_matches, mean(scores), std(scores), ...
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
        sprintf('Prev: %d corners', n_prev), ...
        'Color', [0.7 0.7 1.0], 'FontSize', 9);
    text(hAx, x_offset + 5, nRows + 12, ...
        sprintf('Curr: %d corners', n_curr), ...
        'Color', [0.3 1.0 0.5], 'FontSize', 9);

    hold(hAx, 'off');

    % ----------------------------------------------------------------
    % 7. UI controls
    % ----------------------------------------------------------------
    if ~auto_advance

        % Shared state for button callbacks
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

        % "Stop" button (stop debugging for all remaining frames)
        uicontrol('Parent', hFig, ...
            'Style', 'pushbutton', ...
            'String', 'Stop All', ...
            'FontSize', 10, ...
            'Units', 'normalized', ...
            'Position', [0.56 0.005 0.10 0.055], ...
            'BackgroundColor', [0.5 0.15 0.15], ...
            'ForegroundColor', 'w', ...
            'Callback', @(~,~) onStop(hFig));

        % Frame counter
        uicontrol('Parent', hFig, ...
            'Style', 'text', ...
            'String', sprintf('Frame %d → %d', frameIndex-1, frameIndex), ...
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
%  Helper: map NCC scores to RGB colours
% =====================================================================
function colors = nccToColor(scores, cmap_choice)
%nccToColor Map NCC scores in [0, 1] to RGB colours.
%
%   'rg'     - Red (low NCC) → Yellow (mid) → Green (high)
%   'jet'    - Jet colourmap
%   'parula' - Parula colourmap

    n = numel(scores);
    colors = zeros(n, 3);

    % Clamp to [0, 1]
    s = max(0, min(1, scores));

    switch lower(cmap_choice)
        case 'rg'
            % Red → Yellow → Green
            % R: 1.0 → 1.0 → 0.0   (drops in second half)
            % G: 0.0 → 1.0 → 1.0   (rises in first half)
            % B: 0.0 → 0.0 → 0.0
            for i = 1:n
                if s(i) < 0.5
                    t = s(i) * 2;           % 0 → 1 over first half
                    colors(i,:) = [1.0, t, 0.0];
                else
                    t = (s(i) - 0.5) * 2;   % 0 → 1 over second half
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
            % Fallback to single colour (cyan)
            colors = repmat([0 0.8 0.8], n, 1);
    end
end
