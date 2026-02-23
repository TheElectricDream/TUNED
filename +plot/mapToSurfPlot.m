function [] = mapToSurfPlot(map)

    % Generate Point Vectors from Map
    [pointCloud] = process.generateMeshFromFrame(map);
    
    % Trim the Vectors to Remove NaNs
    x_trimmed = pointCloud(~isnan(pointCloud(:,3)),1);
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)),2);
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)),3);
    
    % Define a Query Grid
    nx = 1000; ny = 1000;
    [xq,yq] = meshgrid(linspace(min(x_trimmed),max(x_trimmed),nx), ...
    linspace(min(y_trimmed),max(y_trimmed),ny));

    % Interpolate Points
    F = scatteredInterpolant(x_trimmed,y_trimmed,...
    z_trimmed,"natural","nearest");
    zq = F(xq,yq);
    
    % Plot the Surface
    figure()
    surf(xq,yq,zq,'EdgeColor','none'); 
    shading interp; 
    view(3);

end
