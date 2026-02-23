function [S, normalized_output_frame] = speedInvariantTimeSurface(S, x_list, y_list, R)
% SPEEDINVARIANTTIMESURFACE Update surface using the paper's specific algorithm.
%
%   [S, normalized_output_frame] = speedInvariantTimeSurface(S, x_list, y_list, R)
%
%   Inputs:
%     S      - The current state of the Surface (matrix of integers)
%     x_list - Vector of row coordinates for incoming events
%     y_list - Vector of column coordinates for incoming events
%     R      - Radius of the neighborhood (scalar integer)
%
%   Outputs:
%     S      - The current state of the Surface (matrix of integers)
%     normalized_output_frame - The normalized surface for visualization
%
%   Notes:
%     - Source paper: https://arxiv.org/pdf/1903.11332

    % Define the maximum value constant
    max_val = (2*R + 1)^2;
    
    % Get image dimensions
    [H, W] = size(S);

    % Iterate through all events in the provided lists
    for k = 1:length(x_list)
        x = x_list(k);
        y = y_list(k);

        % Define the neighborhood bounds (clamping to image edges)
        x_min = max(1, x - R);
        x_max = min(H, x + R);
        y_min = max(1, y - R);
        y_max = min(W, y + R);

        % Get the current value at the event location
        current_val_at_center = S(x, y);

        % Extract the neighborhood
        patch = S(x_min:x_max, y_min:y_max);

        % Find neighbors that need decrementing. 
        % Condition: S(neighbor) >= S(center)
        mask = patch >= current_val_at_center;

        % Decrement those neighbors
        patch(mask) = patch(mask) - 1;

        % Write the patch back
        S(x_min:x_max, y_min:y_max) = patch;

        % Set the center pixel to the maximum value
        S(x, y) = max_val;
    end

    normalized_output_frame = S / max_val;
end