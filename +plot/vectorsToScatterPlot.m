function [] = vectorsToScatterPlot(x, y, z, holdFig)
% VECTORSTOSCATTERPLOT  3D scatter plot from raw coordinate vectors.
%
%   VECTORSTOSCATTERPLOT(X, Y, Z, HOLDFIG) displays a 3D scatter
%   plot colored by z-value.
%
%   Inputs:
%     x, y, z - [N x 1] Coordinate vectors.
%     holdFig - Logical. If true, overlay on existing figure 16.
%
%   See also: plot.mapToScatterPlot

    if holdFig
        figure(16);
        hold on;
    else
        figure();
    end

    scatter3(x, y, z, 36, z, 'filled');
    colormap(jet);
    colorbar;
    view(3);
    set(gca, 'FontSize', 16);

end