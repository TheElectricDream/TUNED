function [] = vectorsToPCViewer(x, y, z)
% VECTORSTOPCVIEWER  Point cloud viewer from raw coordinate vectors.
%
%   VECTORSTOPCVIEWER(X, Y, Z) normalizes the input coordinates
%   and displays them using MATLAB's pcviewer.
%
%   Inputs:
%     x, y, z - [N x 1] Coordinate vectors.
%
%   Notes:
%     - Coordinates are normalized to [0, 1] using the maximum
%       values in each dimension.
%
%   See also: plot.mapToPCViewer

    x_norm = x ./ max(x(:));
    y_norm = y ./ max(y(:));
    z_norm = z ./ max(z(:));

    ptCloudData = pointCloud([x_norm, y_norm, z_norm]);
    pcviewer(ptCloudData);

end