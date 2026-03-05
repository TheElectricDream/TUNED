function [] = mapToScaledImage(map)
% MAPTOSCALEDIMAGE  Display a 2D map with imagesc.
%
%   MAPTOSCALEDIMAGE(MAP) displays the map as a scaled image with
%   axis labels and a colorbar.
%
%   Inputs:
%     map - [H x W] 2D value map.
%
%   See also: imagesc

    figure();
    imagesc(map');
    xlabel('X [px]');
    ylabel('Y [px]');
    colorbar;
    set(gca, 'FontSize', 16);

end