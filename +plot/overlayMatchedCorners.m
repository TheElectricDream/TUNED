function display_frame = overlayMatchedCorners(display_frame, ...
    corners_curr, corners_prev, matches, scores, varargin)
%overlayMatchedCorners Draw match lines between corresponding corners.
%
%   DISPLAY_FRAME = overlayMatchedCorners(DISPLAY_FRAME, CORNERS_CURR,
%       CORNERS_PREV, MATCHES, SCORES)
%   DISPLAY_FRAME = overlayMatchedCorners(..., OPTS)
%
%   Burns lines connecting matched corner pairs and markers at each
%   endpoint directly into a grayscale uint8 display frame.
%
%   Inputs:
%     display_frame - [H x W] uint8 grayscale image (transposed sensor,
%                     as used by the ATS video writer).
%     corners_curr  - [M x 2] current-frame corners [row, col] in
%                     SENSOR coordinates.
%     corners_prev  - [N x 2] previous-frame corners [row, col] in
%                     SENSOR coordinates.
%     matches       - [P x 2] index pairs from matchCVDescriptors:
%                     matches(k,1) → corners_curr index,
%                     matches(k,2) → corners_prev index.
%     scores        - [P x 1] NCC scores per match (used for brightness
%                     scaling: higher NCC → brighter line).
%     opts          - (Optional) struct with fields:
%                       .marker_value - uint8 brightness for endpoints
%                                       (default: 255)
%                       .line_value   - uint8 brightness for lines
%                                       (default: 200)
%                       .arm_length   - crosshair arm [px] (default: 2)
%                       .scale_by_ncc - logical: scale line brightness
%                                       by NCC score (default: true)
%
%   Outputs:
%     display_frame - Modified image with match lines burned in.
%
%   Note on coordinates:
%     Sensor (row, col) → display (col, row) via the transpose.
%
%   See also: features.matchCVDescriptors, plot.overlayHarrisCVMarkers

    % Defaults
    marker_val    = uint8(255);
    line_base_val = 200;
    arm           = 2;
    scale_ncc     = true;

    if nargin >= 6 && isstruct(varargin{1})
        opts = varargin{1};
        if isfield(opts, 'marker_value'), marker_val    = uint8(opts.marker_value); end
        if isfield(opts, 'line_value'),   line_base_val = opts.line_value; end
        if isfield(opts, 'arm_length'),   arm           = opts.arm_length; end
        if isfield(opts, 'scale_by_ncc'), scale_ncc     = opts.scale_by_ncc; end
    end

    if isempty(matches)
        return;
    end

    [H, W] = size(display_frame);

    for k = 1:size(matches, 1)
        % Current corner in display coords
        r_c = corners_curr(matches(k,1), 2);   % sensor col → display row
        c_c = corners_curr(matches(k,1), 1);   % sensor row → display col

        % Previous corner in display coords
        r_p = corners_prev(matches(k,2), 2);
        c_p = corners_prev(matches(k,2), 1);

        % Line brightness (optionally scaled by NCC)
        if scale_ncc
            lv = uint8(round(line_base_val * max(0, scores(k))));
        else
            lv = uint8(line_base_val);
        end

        % Draw line via simple DDA
        display_frame = drawLineDDA(display_frame, ...
            r_c, c_c, r_p, c_p, lv, H, W);

        % Draw marker at current corner
        display_frame = drawCrosshair(display_frame, ...
            r_c, c_c, arm, marker_val, H, W);

        % Draw smaller marker at previous corner (dimmer)
        display_frame = drawCrosshair(display_frame, ...
            r_p, c_p, max(1, arm-1), uint8(round(double(marker_val)*0.6)), H, W);
    end
end


% =====================================================================
%  Local helper: DDA line rasterisation
% =====================================================================
function img = drawLineDDA(img, r1, c1, r2, c2, val, H, W)
    dr = r2 - r1;
    dc = c2 - c1;
    steps = max(abs(dr), abs(dc));
    if steps == 0
        if r1 >= 1 && r1 <= H && c1 >= 1 && c1 <= W
            img(r1, c1) = val;
        end
        return;
    end
    r_inc = dr / steps;
    c_inc = dc / steps;
    r = double(r1);
    c = double(c1);
    for s = 0:steps
        ri = round(r);
        ci = round(c);
        if ri >= 1 && ri <= H && ci >= 1 && ci <= W
            img(ri, ci) = max(img(ri, ci), val);
        end
        r = r + r_inc;
        c = c + c_inc;
    end
end


% =====================================================================
%  Local helper: crosshair marker
% =====================================================================
function img = drawCrosshair(img, r, c, arm, val, H, W)
    c_lo = max(1, c - arm);
    c_hi = min(W, c + arm);
    if r >= 1 && r <= H
        img(r, c_lo:c_hi) = val;
    end
    r_lo = max(1, r - arm);
    r_hi = min(H, r + arm);
    if c >= 1 && c <= W
        img(r_lo:r_hi, c) = val;
    end
end
