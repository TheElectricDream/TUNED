function printPercentComplete(index, total, duration)
% PRINTPERCENTCOMPLETE  Display processing progress to console.
%
%   PRINTPERCENTCOMPLETE(INDEX, TOTAL, DURATION) prints a one-line
%   progress message showing percentage complete and elapsed time.
%
%   Inputs:
%     index    - Current iteration index.
%     total    - Total number of iterations.
%     duration - Elapsed time for the current iteration [s].

    percent_complete = index / total * 100;
    fprintf('Processing... %.1f%% complete (%.2f seconds)\n', ...
        percent_complete, duration);

end