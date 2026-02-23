function [sorted_x, sorted_y, sorted_t, unique_idx, pos, group_ends] = ...
    spreadEventsSpatially(x, y, t, imgSz, radius)
% SPREADEVENTSSPATIALLY Propagates events to neighbors to create a "splatted"
% set of vectors compatible with computeNeighborhoodStats.
%
%   This function effectively performs a spatial convolution on the event
%   stream. By replicating events to neighbors, the 'mean' stats become
%   spatial averages, and 'diff' stats become local coherence measures.
%
% INPUTS:
%   x, y   - Original coordinate vectors (single or double)
%   t      - Original timestamp vector
%   imgSz  - [height, width]
%   radius - The expansion radius (e.g., 2 means a 5x5 block: -2 to +2)
%
% OUTPUTS:
%   All outputs are formatted to be passed directly into 
%   computeNeighborhoodStats.

    % 1. DEFINE KERNEL (Offsets)
    % We create a grid of offsets relative to (0,0)
    % To "Splat" over X by X, radius should be (X-1)/2.
    % E.g. Radius 2 = 5x5 kernel.
    range = -radius:radius;
    [dY, dX] = meshgrid(range, range);
    
    % Flatten offsets
    off_x = int32(dX(:));
    off_y = int32(dY(:));
    num_replicas = length(off_x);

    % 2. REPLICATE DATA (Vectorized)
    % We use indexing to replicate the original vectors K times.
    % This is much faster than loops.
    n_events = length(x);
    
    % Replicate indices: [1,1,1..., 2,2,2..., etc] to keep events grouped first
    % or [1,2,3..., 1,2,3...] to group by offset. 
    % We repeat inner elements to keep memory operations linear.
    % Method: We essentially expand the vectors by factor of K
    
    % Convert to compatible types if needed
    if ~isa(x, 'int32'), x = int32(x); end
    if ~isa(y, 'int32'), y = int32(y); end
    
    % Indexing trick for replication (K times each element)
    % idx_map = repelem(1:n_events, num_replicas)'; % Requires newer MATLAB
    % Faster compatible version:
    idx_map = floor((0:n_events*num_replicas - 1) / num_replicas)' + 1;
    
    % Create the Big Vectors
    x_big = x(idx_map);
    y_big = y(idx_map);
    t_big = t(idx_map);
    
    % Create the Big Offsets
    % We repeat the offset pattern N times
    off_x_big = repmat(off_x, n_events, 1);
    off_y_big = repmat(off_y, n_events, 1);
    
    % Apply offsets
    x_new = x_big + off_x_big;
    y_new = y_big + off_y_big;
    
    % 3. FILTER BOUNDS
    % Remove events pushed off the image
    valid_mask = x_new >= 1 & x_new <= imgSz(1) & ...
                 y_new >= 1 & y_new <= imgSz(2);
             
    x_valid = x_new(valid_mask);
    y_valid = y_new(valid_mask);
    t_valid = t_big(valid_mask);
    
    % 4. SORTING AND GROUPING
    % We must sort by Linear Index (pixel location) AND Timestamp.
    % Sorting by timestamp is critical for 'mean_diff' to work correctly.
    
    linear_idx = sub2ind(imgSz, x_valid, y_valid);
    
    % sortrows is efficient for multi-column sort
    % Col 1: Pixel ID, Col 2: Time
    [~, sort_order] = sortrows([double(linear_idx), double(t_valid)]);
    
    % Apply sort
    sorted_x = x_valid(sort_order);
    sorted_y = y_valid(sort_order);
    sorted_t = t_valid(sort_order);
    sorted_idx = linear_idx(sort_order);
    
    % Generate grouping indices for the stats function
    [unique_idx, pos, ~] = unique(sorted_idx);
    group_ends = [pos(2:end)-1; length(sorted_idx)];

end