function [] = mapToSurfPlot(map)
% MAPTOSURFPLOT  Interpolated surface plot from a 2D value map.
%
%   MAPTOSURFPLOT(MAP) converts a 2D map to a point cloud,
%   interpolates onto a regular grid, and displays as a surface.
%
%   Inputs:
%     map - [H x W] 2D value map.
%
%   See also: process.generateMeshFromFrame, plot.mapToScatterPlot

    [pointCloud] = process.generateMeshFromFrame(map);

    x_trimmed = pointCloud(~isnan(pointCloud(:,3)), 1);
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)), 2);
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)), 3);

    nx = 1000;
    ny = 1000;
    [xq, yq] = meshgrid(...
        linspace(min(x_trimmed), max(x_trimmed), nx), ...
        linspace(min(y_trimmed), max(y_trimmed), ny));

    F = scatteredInterpolant(x_trimmed, y_trimmed, ...
        z_trimmed, 'natural', 'nearest');
    zq = F(xq, yq);

    figure();
    surf(xq, yq, zq, 'EdgeColor', 'none');
    shading interp;
    view(3);
    set(gca, 'FontSize', 16);

end