function [] = vectorsToScatterPlot(x, y, z, holdFig)
    
    % Plot the Surface
    if holdFig
        figure(16)
        scatter3(x, y, z, 36, z, 'filled') 
        hold on;
        colormap(jet); % Add a colormap to the scatter plot
        colorbar; % Optional: Add a colorbar to indicate the color scale
        view(3);
        set(gca, 'FontSize', 16)

    else
        figure()
        scatter3(x, y, z, 36, z, 'filled') 
        colormap(jet); % Add a colormap to the scatter plot
        colorbar; % Optional: Add a colorbar to indicate the color scale
        view(3);
        set(gca, 'FontSize', 16)
    end

end