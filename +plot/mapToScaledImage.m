function [] = mapToScaledImage(map)


    figure();
    imagesc(map');
    xlabel('X [px]');
    ylabel('Y [px]');
    colorbar;
    set(gca, 'FontSize', 16)


end