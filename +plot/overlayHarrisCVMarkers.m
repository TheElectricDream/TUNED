function display_frame = overlayHarrisCVMarkers(display_frame, corners_sub, varargin)
%overlayHarrisCVMarkers Burn Harris-CV corner markers onto a uint8 frame.
%
%   DISPLAY_FRAME = overlayHarrisCVMarkers(DISPLAY_FRAME, CORNERS_SUB)
%   DISPLAY_FRAME = overlayHarrisCVMarkers(DISPLAY_FRAME, CORNERS_SUB, OPTS)
%
%   Burns cross-hair markers at detected Harris-CV corner locations
%   directly into a grayscale uint8 image. The markers are rendered as
%   bright pixels so they are visible on the dark event-camera output.
%
%   This function operates on the DISPLAY frame (transposed: H x W,
%   where H = imgSz(2) and W = imgSz(1)), matching the convention
%   used by the ATS video writer in main.m.
%
%   Inputs:
%     display_frame - [H x W] uint8 grayscale image (transposed sensor)
%     corners_sub   - [K x 2] array of [row, col] in SENSOR coordinates
%                     (imgSz orientation, from detectHarrisCV output).
%                     The function maps (row, col) -> (col, row) in
%                     display coordinates via the transpose.
%     opts          - (Optional) struct with fields:
%                       .marker_value  - uint8 brightness (default: 255)
%                       .arm_length    - cross-hair arm in pixels (default: 3)
%
%   Outputs:
%     display_frame - Modified image with markers burned in.
%
%   Note on coordinate convention:
%     detectHarrisCV returns corners in SENSOR space: row ∈ [1, imgSz(1)],
%     col ∈ [1, imgSz(2)]. The display frame is the transpose of the
%     sensor frame (i.e., display = sensor'), so sensor (row, col) maps
%     to display (col, row).
%
%   See also: features.detectHarrisCV, plot.overlayCornerMarkers

    % Defaults
    marker_val = uint8(255);
    arm        = 3;

    if nargin >= 3 && isstruct(varargin{1})
        opts = varargin{1};
        if isfield(opts, 'marker_value'), marker_val = uint8(opts.marker_value); end
        if isfield(opts, 'arm_length'),   arm = opts.arm_length; end
    end

    if isempty(corners_sub)
        return;
    end

    [H, W] = size(display_frame);

    for i = 1:size(corners_sub, 1)
        % Sensor coordinates -> display coordinates (transpose)
        dr = corners_sub(i, 2);   % sensor col -> display row
        dc = corners_sub(i, 1);   % sensor row -> display col

        % Horizontal arm
        c_lo = max(1, dc - arm);
        c_hi = min(W, dc + arm);
        if dr >= 1 && dr <= H
            display_frame(dr, c_lo:c_hi) = marker_val;
        end

        % Vertical arm
        r_lo = max(1, dr - arm);
        r_hi = min(H, dr + arm);
        if dc >= 1 && dc <= W
            display_frame(r_lo:r_hi, dc) = marker_val;
        end
    end
end
