function [] = mapToPCViewer(map)
% MAPTOPCVIEWER  Point cloud viewer from a 2D value map.
%
%   MAPTOPCVIEWER(MAP) converts a 2D map to a normalized point
%   cloud and displays it using MATLAB's pcviewer.
%
%   Inputs:
%     map - [H x W] 2D value map.
%
%   Notes:
%     - Coordinates are normalized to [0, 1] using imgSz derived
%       from the map dimensions rather than hardcoded values.
%
%   See also: process.generateMeshFromFrame, plot.mapToScatterPlot

    [pointsFromMap] = process.generateMeshFromFrame(map);

    x_trimmed = pointsFromMap(~isnan(pointsFromMap(:,3)), 1);
    y_trimmed = pointsFromMap(~isnan(pointsFromMap(:,3)), 2);
    z_trimmed = pointsFromMap(~isnan(pointsFromMap(:,3)), 3);

    [H, W] = size(map);
    x_norm = x_trimmed ./ W;
    y_norm = y_trimmed ./ H;
    z_norm = z_trimmed ./ max(z_trimmed(:));

    ptCloudData = pointCloud([x_norm, y_norm, z_norm]);
    pcviewer(ptCloudData);

end