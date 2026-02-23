function [] = mapToScatterPlot(map, holdFig)
    % mapToScatterPlot Visualizes a 3D scatter plot from a given map.
    % This function processes the input map to generate a point cloud,
    % trims the vectors to remove NaN values, and then creates a 3D scatter
    % plot with a color gradient based on the z-values.

    % Generate Point Vectors from Map
    [pointCloud] = process.generateMeshFromFrame(map);
    
    % Trim the Vectors to Remove NaNs
    x_trimmed = pointCloud(~isnan(pointCloud(:,3)),1);
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)),2);
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)),3);
    
    % Plot the Surface
    if holdFig
        figure(15)
        scatter3(x_trimmed, y_trimmed, z_trimmed, 36, z_trimmed, 'filled') 
        hold on;
        colormap(jet); % Add a colormap to the scatter plot
        colorbar; % Optional: Add a colorbar to indicate the color scale
        view(3);
        set(gca, 'FontSize', 16)

    else
        figure()
        scatter3(x_trimmed, y_trimmed, z_trimmed, 36, z_trimmed, 'filled') 
        colormap(jet); % Add a colormap to the scatter plot
        colorbar; % Optional: Add a colorbar to indicate the color scale
        view(3);
        set(gca, 'FontSize', 16)
    end

end