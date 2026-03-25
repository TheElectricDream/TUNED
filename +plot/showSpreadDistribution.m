function fig = showSpreadDistribution(map, varargin)
% SHOWSPREADDISTRIBUTION  Pairwise distance histogram for active pixels.
%
%   FIG = SHOWSPREADDISTRIBUTION(MAP) extracts all active pixel
%   coordinates from MAP (map > 0), computes their pairwise Euclidean
%   distances, and displays a histogram of the resulting distribution.
%
%   FIG = SHOWSPREADDISTRIBUTION(MAP, 'Name', Value) accepts optional
%   name-value arguments:
%     'NumBins'     - (100) Number of histogram bins.
%     'FigPosition' - ([100 100 560 420]) Figure position vector.
%
%   Inputs:
%     map - [H x W] 2D value map. Active pixels are those with map > 0.
%
%   Outputs:
%     fig - Figure handle.
%
%   Notes:
%     - Pairwise distance computation via PDIST is O(N^2) in memory and
%       feasible up to approximately 5 000 active pixels. For denser
%       maps, consider subsampling before calling this function.
%     - Coordinates are expressed in pixel units as [col, row].
%
%   See also: plot.mapToPCViewer, plot.mapToScatterPlot

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'NumBins',      100,               ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'FigPosition',  [100 100 560 420],  ...
        @(x) isnumeric(x) && numel(x) == 4);
    parse(p, varargin{:});
    n_bins  = p.Results.NumBins;
    fig_pos = p.Results.FigPosition;

    % ----------------------------------------------------------------
    % 1. Extract active pixel coordinates [col, row]
    % ----------------------------------------------------------------
    [rows, cols] = find(map > 0);
    pts = [cols, rows];

    % ----------------------------------------------------------------
    % 2. Compute pairwise Euclidean distances
    % ----------------------------------------------------------------
    D = pdist(pts, 'euclidean');   % upper-triangle vector, length N*(N-1)/2

    % ----------------------------------------------------------------
    % 3. Plot histogram
    % ----------------------------------------------------------------
    fig = figure('Position', fig_pos);
    histogram(D, n_bins);
    xlabel('Distance (px)');
    ylabel('Count');
    title('Pairwise Spread Distribution');
    grid on;

end