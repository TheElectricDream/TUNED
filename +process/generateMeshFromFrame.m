function [pointCloud] = generateMeshFromFrame(map)
% GENERATEMESHFROMFRAME  Convert 2D map to Nx3 point cloud.
%
%   POINTCLOUD = GENERATEMESHFROMFRAME(MAP) creates a point cloud
%   from a 2D value map, where x and y are pixel coordinates and
%   z is the map value. Zero-valued pixels are set to NaN.
%
%   Inputs:
%     map - [H x W] 2D matrix of values (e.g., time surface).
%
%   Outputs:
%     pointCloud - [H*W x 3] Array of [x, y, z] coordinates.
%                  NaN rows correspond to zero-valued pixels.
%
%   See also: plot.mapToScatterPlot, plot.mapToSurfPlot

    map(map == 0) = nan;

    [x, y] = meshgrid(1:size(map, 2), 1:size(map, 1));

    pointCloud = [x(:), y(:), map(:)];

end