function [S, normalized_output_frame] = ...
    speedInvariantTimeSurface(S, x_list, y_list, R)
% SPEEDINVARIANTTIMESURFACE  Speed Invariant Time Surface (SITS).
%
%   [S, NORMALIZED_OUTPUT_FRAME] = SPEEDINVARIANTTIMESURFACE(S,
%   X_LIST, Y_LIST, R) updates the surface using the event-driven
%   decrement strategy. Each incoming event sets its pixel to the
%   maximum value and decrements neighbors that have equal or higher
%   values. This produces an edge gradient whose slope is constant
%   regardless of edge speed.
%
%   Reference:
%       Manderscheid et al., "Speed Invariant Time Surface for
%       Learning to Detect Corner Points with Event-Based Cameras,"
%       Proc. IEEE CVPR, pp. 10245-10254, 2019.
%       DOI: 10.1109/CVPR.2019.01049
%
%   Inputs:
%     S      - [imgSz] Current surface state (integer-valued).
%              Initialize with zeros(imgSz).
%     x_list - [N x 1] Row coordinates of events.
%     y_list - [N x 1] Column coordinates of events.
%     R      - Scalar integer neighborhood half-size [pixels].
%              Observation window is (2R+1) x (2R+1).
%
%   Outputs:
%     S                       - [imgSz] Updated surface (integers).
%     normalized_output_frame - [imgSz] Display surface in [0, 1].
%
%   Algorithm:
%     For each event at pixel (x, y):
%       1. Extract the (2R+1)x(2R+1) neighborhood patch.
%       2. Decrement all neighbors with value >= S(x, y).
%       3. Set S(x, y) to the maximum value (2R+1)^2.
%
%   Notes:
%     - Event-driven decrement (not time-driven decay). Values are
%       integer counts, not exponential recency.
%     - Events must be processed sequentially (order matters).
%     - Unsigned output in [0, 1]. Normalization: S / max_val.
%     - Known limitation: lingering problem when motion direction
%       changes (see Xu et al., IJCV 2025, Fig. 4).
%     - Coordinates: x = row, y = col.
%
%   See also: accumulator.timeSurface,
%             accumulator.motionEncodedTimeSurface

    % ----------------------------------------------------------------
    % 0. Constants
    % ----------------------------------------------------------------
    MAX_VAL = (2*R + 1)^2;
    [H, W] = size(S);

    % ----------------------------------------------------------------
    % 1. Process events sequentially (order-dependent)
    % ----------------------------------------------------------------
    for k = 1:length(x_list)
        x = x_list(k);
        y = y_list(k);

        % Clamp neighborhood to image bounds
        x_min = max(1, x - R);
        x_max = min(H, x + R);
        y_min = max(1, y - R);
        y_max = min(W, y + R);

        % Extract neighborhood patch
        current_val = S(x, y);
        patch = S(x_min:x_max, y_min:y_max);

        % Decrement neighbors with value >= center
        mask = patch >= current_val;
        patch(mask) = patch(mask) - 1;

        % Write patch back and set center to max
        S(x_min:x_max, y_min:y_max) = patch;
        S(x, y) = MAX_VAL;
    end

    % ----------------------------------------------------------------
    % 2. Normalize for display
    % ----------------------------------------------------------------
    % Count-based surface [0, max_val] → [0, 1].
    normalized_output_frame = S / MAX_VAL;

end