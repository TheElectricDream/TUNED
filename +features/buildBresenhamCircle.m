function offsets = buildBresenhamCircle(radius)
%buildBresenhamCircle Generate ordered pixel offsets for a Bresenham circle.
%
%   OFFSETS = buildBresenhamCircle(RADIUS) returns an N×2 matrix of [dx, dy]
%   offsets arranged in clockwise order starting from (0, -radius), where
%   the offsets describe a discrete Bresenham midpoint circle of the given
%   RADIUS.
%
%   This is used to construct the circular inspection masks required by the
%   Arc* corner detector (Alzugaray & Chli, IEEE RA-L 2018). The ordering
%   is consistent with the FAST-style convention used in eFAST and Arc*.
%
%   Inputs:
%       radius  - Positive integer circle radius (typically 3 or 4)
%
%   Outputs:
%       offsets - [N × 2] matrix of [dx, dy] pixel offsets, ordered CW
%
%   Reference:
%       Alzugaray, I. and Chli, M. (2018), "Asynchronous Corner Detection
%       and Tracking for Event Cameras in Real-Time," IEEE Robotics and
%       Automation Letters, 3(4), pp. 3177–3184.
%
%   See also: features.detectArcStarCorners

    % --- Midpoint circle algorithm (first octant) ---
    % Generate all unique points on the discrete circle, then sort CW.
    x = 0;
    y = radius;
    d = 1 - radius;  % Decision parameter

    pts = [];  % Collector for octant-reflected points

    while x <= y
        % Eight-way symmetry reflections
        pts = [pts;
            x,  y;
           -x,  y;
            x, -y;
           -x, -y;
            y,  x;
           -y,  x;
            y, -x;
           -y, -x]; %#ok<AGROW>

        x = x + 1;
        if d < 0
            d = d + 2*x + 1;
        else
            y = y - 1;
            d = d + 2*(x - y) + 1;
        end
    end

    % Remove duplicate points (octant boundaries produce overlaps)
    pts = unique(pts, 'rows', 'stable');

    % --- Sort clockwise from (0, -radius) using atan2 ---
    % atan2 convention: angle measured from +x axis, CCW positive.
    % We want CW from "north" = (0, -radius).
    % Map: theta_cw = atan2(dx, -dy)  [CW angle from north]
    angles = atan2(pts(:,1), -pts(:,2));

    % Wrap to [0, 2*pi) for a clean CW sort starting from north
    angles = mod(angles, 2*pi);

    [~, sort_idx] = sort(angles);
    offsets = pts(sort_idx, :);

end
