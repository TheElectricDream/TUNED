function rgb_frame = overlayCornerMarkers(gray_frame, corner_events, marker_opts)
%overlayCornerMarkers Draw corner-event markers on a grayscale frame.
%
%   RGB_FRAME = overlayCornerMarkers(GRAY_FRAME, CORNER_EVENTS) converts
%   the grayscale frame to RGB and stamps a small colored cross at each
%   detected corner location. Positive-polarity corners are drawn in green,
%   negative-polarity corners in magenta.
%
%   RGB_FRAME = overlayCornerMarkers(GRAY_FRAME, CORNER_EVENTS, MARKER_OPTS)
%   uses custom marker settings from the MARKER_OPTS struct.
%
%   This function expects the frame to ALREADY BE TRANSPOSED into display
%   orientation (H x W), matching the convention in main.m where frames
%   are transposed before being passed to writeVideo and imwrite:
%
%       grayscale_normalized_output_frame'   % <-- H x W display orientation
%
%   Corner event coordinates (x, y) from detectArcStarCorners use the
%   IEI-ATS convention where x = row and y = column in imgSz = [nRows,
%   nCols] space. After the transpose, these map to display pixel (y, x),
%   i.e. column y becomes the display row, and row x becomes the display
%   column. This mapping is handled internally.
%
%   Inputs:
%       gray_frame    - [H x W] uint8 grayscale frame (display orientation,
%                       i.e. already transposed from imgSz layout)
%       corner_events - [K x 4] array of corners [x, y, t, polarity] from
%                       features.detectArcStarCorners. May be empty.
%       marker_opts   - (Optional) struct with fields:
%                       .radius     - Cross arm length in pixels (default 2)
%                       .color_pos  - [R G B] uint8 for positive polarity
%                                     (default [0 255 0] — green)
%                       .color_neg  - [R G B] uint8 for negative polarity
%                                     (default [255 0 255] — magenta)
%                       .thickness  - Line thickness in pixels (default 1)
%
%   Output:
%       rgb_frame - [H x W x 3] uint8 RGB frame with markers overlaid
%
%   See also: features.detectArcStarCorners

    % --- Default marker options ---
    if nargin < 3 || isempty(marker_opts)
        marker_opts = struct();
    end
    if ~isfield(marker_opts, 'radius'),    marker_opts.radius    = 2;             end
    if ~isfield(marker_opts, 'color_pos'), marker_opts.color_pos = uint8([0 255 0]);   end
    if ~isfield(marker_opts, 'color_neg'), marker_opts.color_neg = uint8([255 0 255]); end
    if ~isfield(marker_opts, 'thickness'), marker_opts.thickness = 1;             end

    radius    = marker_opts.radius;
    color_pos = marker_opts.color_pos;
    color_neg = marker_opts.color_neg;
    thick     = marker_opts.thickness;

    % --- Promote grayscale to RGB ---
    [H, W] = size(gray_frame);
    rgb_frame = repmat(gray_frame, [1, 1, 3]);

    % --- Early exit if no corners ---
    if isempty(corner_events) || size(corner_events, 1) == 0
        return;
    end

    % --- Extract corner coordinates and polarity ---
    % Columns: [x_row, y_col, t, polarity]
    cx_row = corner_events(:, 1);  % Row index in imgSz space
    cy_col = corner_events(:, 2);  % Col index in imgSz space
    cp     = corner_events(:, 4);  % Polarity (-1 or +1)

    % Map from imgSz [nRows x nCols] to transposed display [nCols x nRows]:
    %   display_row = cy_col    (original column becomes display row)
    %   display_col = cx_row    (original row becomes display column)
    dr = cy_col;   % display row
    dc = cx_row;   % display col

    % --- Build cross marker offset pattern ---
    % A cross with arms of length 'radius' and width 'thickness'
    offsets = [];
    half_t = floor((thick - 1) / 2);

    % Horizontal arm
    for dy = -radius:radius
        for dx = -half_t:half_t
            offsets = [offsets; dx, dy]; %#ok<AGROW>
        end
    end

    % Vertical arm
    for dx = -radius:radius
        for dy = -half_t:half_t
            offsets = [offsets; dx, dy]; %#ok<AGROW>
        end
    end

    % Remove duplicates (center pixel overlap)
    offsets = unique(offsets, 'rows');

    n_offsets = size(offsets, 1);
    n_corners = numel(dr);

    % --- Stamp markers ---
    for k = 1:n_corners
        % Select color based on polarity
        if cp(k) > 0
            clr = color_pos;
        else
            clr = color_neg;
        end

        for j = 1:n_offsets
            pr = dr(k) + offsets(j, 1);
            pc = dc(k) + offsets(j, 2);

            % Bounds check
            if pr >= 1 && pr <= H && pc >= 1 && pc <= W
                rgb_frame(pr, pc, 1) = clr(1);
                rgb_frame(pr, pc, 2) = clr(2);
                rgb_frame(pr, pc, 3) = clr(3);
            end
        end
    end

end