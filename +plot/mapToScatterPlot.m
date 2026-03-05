function [] = mapToScatterPlot(map, holdFig)
% MAPTOSCATTERPLOT  3D scatter plot from a 2D value map.
%
%   MAPTOSCATTERPLOT(MAP, HOLDFIG) converts a 2D map to a point
%   cloud and displays it as a 3D scatter plot colored by z-value.
%
%   Inputs:
%     map     - [H x W] 2D value map.
%     holdFig - Logical. If true, overlay on existing figure 15.
%
%   See also: process.generateMeshFromFrame, plot.mapToSurfPlot

    [pointCloud] = process.generateMeshFromFrame(map);

    x_trimmed = pointCloud(~isnan(pointCloud(:,3)), 1);
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)), 2);
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)), 3);

    if holdFig
        figure(15);
        hold on;
    else
        figure();
    end

    scatter3(x_trimmed, y_trimmed, z_trimmed, ...
        36, z_trimmed, 'filled');
    colormap(jet);
    colorbar;
    view(3);
    set(gca, 'FontSize', 16);

end