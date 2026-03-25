function [] = visualizePersistance(t_row, t_col, t_val, ...
    b_rows, b_cols, b_vals, sortedIdx)
% VISUALIZEPERSISTANCE  3D plot of persistence neighbourhood search.
%
%   VISUALIZEPERSISTANCE(T_ROW, T_COL, T_VAL, B_ROWS, B_COLS,
%   B_VALS, SORTEDIDX) visualizes a target point from mapA and its
%   nearest neighbours from mapB, with lines connecting them.
%
%   Inputs:
%     t_row, t_col, t_val - Scalar target point coordinates and
%                           value (from mapA).
%     b_rows, b_cols      - [M x 1] Row/col indices of mapB points.
%     b_vals              - [M x 1] Values of mapB points.
%     sortedIdx           - [K x 1] Indices into b_rows/b_cols of
%                           the closest neighbours.
%
%   Notes:
%     - Creates a new figure. Context points are blue, closest
%       points are green, target is a red star.
%
%   See also: coherence.findPersistenceVectorized

    figure('Color', 'white');
    clf; hold on; grid on; rotate3d on;

    % Separate closest from rest
    isClosest = false(size(b_rows));
    isClosest(sortedIdx) = true;

    rest_rows = b_rows(~isClosest);
    rest_cols = b_cols(~isClosest);
    rest_vals = b_vals(~isClosest);

    closest_rows = b_rows(isClosest);
    closest_cols = b_cols(isClosest);
    closest_vals = b_vals(isClosest);

    % Context points (small, transparent blue)
    active = rest_vals ~= 0;
    scatter3(rest_cols(active), rest_rows(active), ...
        rest_vals(active), 10, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.2);

    % Closest points (large green circles)
    scatter3(closest_cols, closest_rows, closest_vals, ...
        100, 'g', 'filled', 'MarkerEdgeColor', 'k');

    % Target point (red star)
    scatter3(t_col, t_row, t_val, 200, 'r', 'p', ...
        'filled', 'MarkerEdgeColor', 'k');

    % Connecting lines
    for i = 1:length(closest_rows)
        plot3([t_col, closest_cols(i)], ...
            [t_row, closest_rows(i)], ...
            [t_val, closest_vals(i)], ...
            'k-', 'LineWidth', 0.5);
    end

    xlabel('Column Index (X)');
    ylabel('Row Index (Y)');
    zlabel('Time Value');
    title('3D Proximity Search: Space + Time');
    view(3);
    set(gca, 'FontSize', 16);
    hold off;

end