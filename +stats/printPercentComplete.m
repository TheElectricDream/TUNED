function printPercentComplete(index, total, duration_filter, duration_accumulator)
% PRINTPERCENTCOMPLETE  Display processing progress to console.
%
%   PRINTPERCENTCOMPLETE(INDEX, TOTAL, DURATION_FILTER, DURATION_ACCUMULATOR) 
%   prints a one-line progress message showing percentage complete and 
%   elapsed time.
%
%   Inputs:
%     index    - Current iteration index.
%     total    - Total number of iterations.
%     duration_filter - Elapsed time for the current filtering iteration [s].
%     duration_accumulator - Elapsed time for the current accumulator iteration [s].
    percent_complete = index / total * 100;
    fprintf('\nProcessing... %.1f%% complete. Stats: \n \t Filtering: (%.3f seconds)\n  \t Accumulation: (%.3f seconds)\n', ...
        percent_complete, duration_filter, duration_accumulator);

end