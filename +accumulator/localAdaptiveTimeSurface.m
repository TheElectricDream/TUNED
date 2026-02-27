function [normalized_output_frame, time_surface_map_raw, tau_filtered, adaptive_gains] = ...
    localAdaptiveTimeSurface(t_mean, time_surface_map_prev, alts_params,...
    filter_mask, polarity_map, counts)

    % Extract parameters
    surface_tau_release  = alts_params.surface_tau_release;
    dt                   = alts_params.dt;
    filter_sigma         = alts_params.filter_sigma;
    filter_size          = alts_params.filter_size;
    div_norm_exp         = alts_params.div_norm_exp;
    symmetric_tone_scale = alts_params.symmetric_tone_scale;

    % Set any NaN values to 0 for computation
    t_mean(isnan(t_mean)) = 0; 

    % Calculate the candidate decay time
    tau_active = max(t_mean, eps);
    
    % Smooth out the active tau
    tau_filtered = imgaussfilt(tau_active, filter_sigma,...
        "FilterSize", filter_size); 

    % Apply the filter_mask BEFORE computing activity so the
    % attack-release envelope sees only events that survive filtering.
    % Without this, filtered pixels appear "active" and get the slow
    % tau_active decay instead of the fast tau_release, causing trails.
    masked_input = polarity_map .* filter_mask;

    % Build activity indicator from the FILTERED input
    activity_indicator = single(masked_input ~= 0);
    activity_indicator(isnan(activity_indicator))=0;

    % Smooth it spatially so the transition isn't pixel-sharp
    se = strel('disk', 1); 
    activity_blurred = imdilate(activity_indicator, se);

    % Active pixels keep tau_active, idle pixels with residual get tau_release
    tau_effective = tau_active .* activity_blurred + ...
                    surface_tau_release .* (1 - activity_blurred);

    % Compute per-pixel blending coefficient (coupled gain + decay)
    adaptive_gains = 1 - exp(-dt ./ tau_effective);
    adaptive_gains(activity_blurred == 0) = 0;
    adaptive_gains = imgaussfilt(adaptive_gains, filter_sigma, ...
         "FilterSize", filter_size);

    % EMA update 
    time_surface_map_raw = adaptive_gains .* masked_input ...
                         + (1 - adaptive_gains) .* time_surface_map_prev;
    time_surface_map_raw = time_surface_map_raw.*single(activity_blurred>0);
    
    % Extract the magnitude of the time_surface_map
    magnitude = abs(time_surface_map_raw);

    % Get the activity pool
    activity_pool = imgaussfilt(magnitude, filter_sigma, 'FilterSize',...
        filter_size);

    % Build normalization pool from event counts
    counts_smooth = imgaussfilt(single(counts), filter_sigma, 'FilterSize',...
        filter_size);
    
    % Divisive normalization: sigma controls the crossover, so regions with 
    % counts >> sigma get compressed toward 1/counts_smooth and regions 
    % with counts << sigma pass through nearly unchanged.
    time_surface_map    = zeros(size(masked_input));
    good_mask           = (abs(masked_input)>0);
    sigma               = median(activity_pool(abs(activity_pool)>0));         

    % Divisive normalization: signed surface / unsigned activity
    time_surface_map(good_mask) = time_surface_map_raw(good_mask) ./...
        (sigma + counts_smooth(good_mask) .^ div_norm_exp);
    
    % Simple outlier rejection for some final cleanup
    mean_value_pos = mean(time_surface_map(time_surface_map>0));
    mean_value_neg = mean(time_surface_map(time_surface_map<0));
    std_value_pos = std(time_surface_map(time_surface_map>0));
    std_value_neg = std(time_surface_map(time_surface_map<0));

    % Reject points which are still sigma outside the mean
    % Set a threshold
    pos_threshold = mean_value_pos + 4*std_value_pos;
    neg_threshold = mean_value_neg - 4*std_value_neg;
    
    % Reject the points
    time_surface_map(time_surface_map>pos_threshold) = ...
        median(time_surface_map(time_surface_map>0));
    time_surface_map(time_surface_map<neg_threshold) = ...
        median(time_surface_map(time_surface_map<0));
    time_surface_map_raw(time_surface_map>pos_threshold) = ...
        median(time_surface_map_raw(time_surface_map>0));
    time_surface_map_raw(time_surface_map<neg_threshold) = ...
        median(time_surface_map_raw(time_surface_map<0));

    % Normalize the output frame
    normalized_output_frame = ...
        process.symmetricToneMappingNorm(time_surface_map,...
        symmetric_tone_scale);

end

