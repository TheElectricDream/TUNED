function [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, current_idx)
    
    %SLICETOVALIDRANGE Slices the event data based on time and filters
    % the results to only include valid events.
    %
    % Inputs:
    %   t_range_n   - A scalar representing the upper limit of the time range.
    %   xk          - A vector of x-coordinates.
    %   yk          - A vector of y-coordinates.
    %   tk          - A vector of time values corresponding to the coordinates.
    %   pk          - A vector of the polarity
    %   current_idx - The index from which to start searching for valid data.
    %
    % Outputs:
    %   current_idx - The updated index after processing the data.
    %   x_valid     - A vector of x-coordinates that are within the valid range.
    %   y_valid     - A vector of y-coordinates that are within the valid range.
    %   t_valid     - A vector of time values that correspond to the valid coordinates.
    %   p_valid     - A vector of polarity values that are within the valid
    %   range.

    % Find the slice of data efficiently by only searching from where we
    % ended in the last loop
    start_node = current_idx;
    while current_idx <= length(tk) && tk(current_idx) <= t_range_n
        current_idx = current_idx+1;
    end
    end_node = current_idx-1;
    
    % Slice the data according to the valid range
    if end_node >= start_node
        x_slice = xk(start_node:end_node);
        y_slice = yk(start_node:end_node);
        t_slice = tk(start_node:end_node);
        p_slice = pk(start_node:end_node);

        % Filter data to the valid spacial range
        valid_mask = x_slice >= 1 & x_slice <= imgSz(1) & ...
                     y_slice >= 1 & y_slice <= imgSz(2);
        
        x_valid = x_slice(valid_mask);
        y_valid = y_slice(valid_mask);
        t_valid = t_slice(valid_mask);
        p_valid = p_slice(valid_mask);

    else
        % Mark as empty arrays
        x_valid = [];
        y_valid = [];
        t_valid = [];
        p_valid = [];
    end

end