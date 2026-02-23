function [] = visualizePersistance(t_row, t_col, t_val, b_rows, b_cols, b_vals, sortedIdx)
% Visualize persistence: plots a target point from mapA and its nearest neighbors from mapB
% Usage:
%   visualizePersistance(t_row, t_col, t_val, ...
%       b_rows, b_cols, b_vals, sortedIdx)
%
% Inputs:
%   t_row, t_col, t_val   - scalar row, column, and time/value of the target point (from mapA)
%   b_rows, b_cols        - vectors of row and column indices for points in mapB
%   b_vals                - vector of time/value for points in mapB
%   sortedIdx             - linear indices (into b_rows/b_cols/b_vals) marking the "closest" points
%
% The function creates a 3D scatter showing the context points (mapB), the top-10 closest
% points highlighted, and the target point from mapA. Lines connect the target to each closest point.

% Create a Figure
figure('Color', 'white');
clf; hold on; grid on; rotate3d on;

% Separate "Closest" from "Rest"
% Create a logical mask for the top 10 indices
isClosest = false(size(b_rows));
isClosest(sortedIdx) = true;

% Split the data
rest_rows = b_rows(~isClosest);
rest_cols = b_cols(~isClosest);
rest_vals = b_vals(~isClosest);

closest_rows = b_rows(isClosest);
closest_cols = b_cols(isClosest);
closest_vals = b_vals(isClosest);

% Plot "Rest of mapB" (Context)
% Style: Small blue dots, high transparency to reduce clutter
scatter3(rest_cols(rest_vals~=0), rest_rows(rest_vals~=0), rest_vals(rest_vals~=0), 10, 'b', 'filled', ...
    'MarkerFaceAlpha', 0.2);

% Plot "10 Closest Points"
% Style: Large green circles, solid
scatter3(closest_cols, closest_rows, closest_vals, 100, 'g', 'filled', ...
    'MarkerEdgeColor', 'k');

% Plot "Target Point from mapA"
% Style: Large red pentagram (star)
scatter3(t_col, t_row, t_val, 200, 'r', 'p', 'filled', ...
    'MarkerEdgeColor', 'k');

% Formatting
xlabel('Column Index (X)');
ylabel('Row Index (Y)');
zlabel('Time Value');
title('3D Proximity Search: Space + Time');
view(3);

% Draw lines connecting target to the 10 closest
% This helps visually measure the "distance"
for i = 1:length(closest_rows)
    plot3([t_col, closest_cols(i)], ...
        [t_row, closest_rows(i)], ...
        [t_val, closest_vals(i)], 'k-', 'LineWidth', 0.5);
end
set(gca, 'FontSize', 16);
hold off;

end
