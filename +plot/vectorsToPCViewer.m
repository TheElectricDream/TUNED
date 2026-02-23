function [] = vectorsToPCViewer(x, y, z)
    % mapToScatterPlot Visualizes a 3D scatter plot from a given map.
    % This function processes the input map to generate a point cloud,
    % trims the vectors to remove NaN values, and then creates a 3D scatter
    % plot with a color gradient based on the z-values.

    % Normalize for easier viewing
    x_norm = x./640;
    y_norm = y./480;
    z_norm = z./max(z(:));

    ptCloudData = pointCloud([x_norm, y_norm, z_norm]);
    
    % Plot the point cloud
    pcviewer(ptCloudData);

end