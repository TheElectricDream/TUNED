function [pointCloud] = generateMeshFromFrame(map)
    % generateMeshFromFrame creates a point cloud from a given 2D map.
    % Inputs:
    %   map - A 2D matrix representing the height values (z-coordinates) for the point cloud.
    % Outputs:
    %   pointCloud - A Nx3 matrix where N is the number of points, containing the x, y, and z coordinates.

    % To keep things clean, set all zeros to NaNs
    map(map==0)=nan;

    % Create meshgrid for x and y coordinates
    [x, y] = meshgrid(1:size(map, 2), 1:size(map, 1));
    
    % Extract z values from tk_diff_map
    z = map;
    
    % Create point cloud
    pointCloud = [x(:), y(:), z(:)];
end