function printPercentComplete(index, total, duration)
    % printPercentComplete - Displays the percentage of completion
    % 
    % Syntax: printPercentComplete(index, total)
    %
    % Inputs:
    %    index - Current index of the process
    %    total - Total number of items to process
    %
    % This function calculates the percentage of completion based on the
    % current index and total items, and prints it to the console.

    percent_complete = index / total * 100; % Calculate the percentage complete
    fprintf('Processing... %.1f%% complete (%.2f seconds)\n',...
        percent_complete, duration); % Display the result
end