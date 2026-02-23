function [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
    computeNeighborhoodStats(sorted_t, unique_idx, pos, group_ends, imgSz)
    
    % Initialize
    t_mean = zeros(imgSz);
    t_std  = zeros(imgSz);
    t_max = zeros(imgSz);
    t_min = zeros(imgSz);
    t_mean_diff = zeros(imgSz);
    t_std_diff = zeros(imgSz);
    
    % Extract each chunk sequentially and calculate the
    % stats of that "chunk".
    for k = 1:length(unique_idx)
        val_chunk_t = sorted_t(pos(k):group_ends(k));
        idx = unique_idx(k);
        t_mean(idx) = mean(val_chunk_t);
        t_max(idx)  = max(val_chunk_t);
        t_min(idx)  = min(val_chunk_t);
        t_std(idx)  = std(val_chunk_t);
    
        if numel(val_chunk_t) > 1
            d = diff(val_chunk_t);
            t_mean_diff(idx) = mean(d);
            t_std_diff(idx)  = std(d);
        else
            t_mean_diff(idx) = 0;   % no interval observable
            t_std_diff(idx)  = 0;
        end
    end

end
