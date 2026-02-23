function [] = mapToPCViewer(map)
    % mapToScatterPlot Visualizes a 3D scatter plot from a given map.
    % This function processes the input map to generate a point cloud,
    % trims the vectors to remove NaN values, and then creates a 3D scatter
    % plot with a color gradient based on the z-values.

    % Generate Point Vectors from Map
    [pointsFromMap] = process.generateMeshFromFrame(map);
    
    % Trim the Vectors to Remove NaNs
    x_trimmed = pointsFromMap(~isnan(pointsFromMap(:,3)),1);
    y_trimmed = pointsFromMap(~isnan(pointsFromMap(:,3)),2);
    z_trimmed = pointsFromMap(~isnan(pointsFromMap(:,3)),3);

    % Normalize for easier viewing
    x_norm = x_trimmed./640;
    y_norm = y_trimmed./480;
    z_norm = z_trimmed./max(z_trimmed(:));

    ptCloudData = pointCloud([x_norm, y_norm, z_norm]);
    
    % Plot the point cloud
    pcviewer(ptCloudData);

end