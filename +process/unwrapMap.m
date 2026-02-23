function [filtered_events] = unwrapMap(map, sorted_x, sorted_y, sorted_t, t_offset) %#ok<*DEFNU>
    % unwrapCoherenceMap - Extract events that pass coherence filtering
    %
    % Inputs:
    %   coherence_map - [H x W] coherence score map
    %   sorted_x      - [N x 1] x-coordinates of events in current frame
    %   sorted_y      - [N x 1] y-coordinates of events in current frame  
    %   sorted_t      - [N x 1] timestamps of events (relative to frame start)
    %   t_offset      - scalar offset to convert to absolute time
    %
    % Outputs:
    %   filtered_events - [M x 4] array of [t, x, y, coherence_score]
    
    % Initialize output
    filtered_events = [];
    
    % For each event, check if its pixel location has sufficient coherence
    for i = 1:length(sorted_x)
        x = sorted_x(i);
        y = sorted_y(i);
        t = sorted_t(i);
        
        % Look up coherence score at this pixel
        coherence_score = map(y, x);
        
        % Only keep events that pass threshold (coherence_map already filtered)
        if coherence_score > 0  % Already thresholded in your main code
            filtered_events = [filtered_events; x, y,  t + t_offset, coherence_score]; %#ok<AGROW>
        end
    end
end