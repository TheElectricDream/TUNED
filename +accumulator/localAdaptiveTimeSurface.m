function [normalized_output_frame, time_surface_map_raw, ...
    tau_filtered, adaptive_gains] = ...
    localAdaptiveTimeSurface(t_mean, time_surface_map_prev, ...
    alts_params, filter_mask, polarity_map, counts)
% LOCALADAPTIVETIMESURFACE  IEI-ATS adaptive local time surface.
%
%   [NORMALIZED_OUTPUT_FRAME, TIME_SURFACE_MAP_RAW, TAU_FILTERED,
%   ADAPTIVE_GAINS] = LOCALADAPTIVETIMESURFACE(T_MEAN,
%   TIME_SURFACE_MAP_PREV, ALTS_PARAMS, FILTER_MASK, POLARITY_MAP,
%   COUNTS) computes the Adaptive Local Time Surface (ALTS) using an
%   IIR/EMA accumulator with per-pixel adaptive decay driven by local
%   inter-event-interval statistics, an asymmetric attack-release
%   envelope, and Carandini-Heeger divisive normalization.
%
%   This is the core contribution of the IEI-ATS algorithm.
%
%   Inputs:
%     t_mean               - [imgSz] Per-pixel EMA-smoothed mean IEI.
%                            Used directly as the active decay time
%                            constant. NaN values are treated as zero.
%     time_surface_map_prev - [imgSz] Previous frame's raw (unnorm.)
%                            surface state. Feeds back into the EMA.
%     alts_params          - Struct with fields:
%       .dt                  - Frame interval [s] (numerator for alpha)
%       .filter_sigma        - Gaussian smoothing sigma [pixels]
%       .filter_size         - Gaussian kernel size [pixels, odd]
%       .surface_tau_release - Idle-pixel release time constant [s]
%       .div_norm_exp        - Divisive normalization exponent gamma
%       .symmetric_tone_scale - Sigmoid tone mapping scale parameter
%     filter_mask          - [imgSz] Coherence filter mask (0 or 1).
%                            Events at zero-mask pixels are excluded.
%     polarity_map         - [imgSz] Signed polarity accumulation for
%                            current frame. Typically computed as:
%                            accumarray([x,y], p_signed, imgSz, @sum, 0)
%     counts               - [imgSz] Per-pixel event count for the
%                            current frame (used in divisive norm.).
%
%   Outputs:
%     normalized_output_frame - [imgSz] Display surface in [0, 1],
%                               after divisive norm. + tone mapping.
%     time_surface_map_raw    - [imgSz] Raw bipolar EMA surface
%                               (unnormalized). Feed back as prev.
%     tau_filtered            - [imgSz] Spatially smoothed active tau.
%     adaptive_gains          - [imgSz] Per-pixel blending coefficient
%                               alpha_eff after spatial smoothing.
%
%   Algorithm:
%     1. Map IEI directly to active tau: tau = max(IEI, eps).
%     2. Spatially smooth tau with a Gaussian kernel.
%     3. Build activity indicator from coherence-filtered polarity,
%        expanded via morphological dilation (disk, r=1).
%     4. Blend attack/release: active pixels use tau_active, idle
%        pixels use tau_release.
%     5. Compute coupled gain-decay: alpha = 1 - exp(-dt/tau_eff).
%     6. EMA update: S = alpha*u + (1-alpha)*S_prev.
%     7. Suppress inactive pixels (zero outside activity mask).
%     8. Divisive normalization: S / (sigma^gamma + C^gamma).
%     9. Tone mapping via symmetric sigmoid.
%
%   Notes:
%     - IIR/EMA-based: recursive update, NOT reset-based. The surface
%       preserves signed polarity (+/- events distinguishable).
%     - Bipolar output (pre-normalization) in approximately [-1, +1].
%       Tone-mapped output in [0, 1] with 0.5 = zero.
%     - Raw unnormalized surface must feed back into next frame
%       (not the normalized output — this prevents self-sustaining
%       residuals from divisive normalization feedback).
%     - Coordinates: x = row, y = col.
%
%   See also: accumulator.timeSurface,
%             accumulator.adaptiveGlobalDecay,
%             coherence.computeCoherenceMask,
%             process.symmetricToneMappingNorm

    % ----------------------------------------------------------------
    % 0. Parse parameters
    % ----------------------------------------------------------------
    surface_tau_release  = alts_params.surface_tau_release;
    dt                   = alts_params.dt;
    filter_sigma         = alts_params.filter_sigma;
    filter_size          = alts_params.filter_size;
    div_norm_exp         = alts_params.div_norm_exp;
    symmetric_tone_scale = alts_params.symmetric_tone_scale;

    % Sanitize IEI input
    t_mean(isnan(t_mean)) = 0;

    % ----------------------------------------------------------------
    % 1. Time constant mapping (identity: tau = IEI)
    % ----------------------------------------------------------------
    tau_active = max(t_mean, eps);

    % Spatial smoothing to prevent pixel-sharp transitions
    tau_filtered = imgaussfilt(tau_active, filter_sigma, ...
        'FilterSize', filter_size);

    % ----------------------------------------------------------------
    % 2. Activity indicator (morphological dilation)
    % ----------------------------------------------------------------
    % Apply coherence mask BEFORE computing activity so the envelope
    % sees only events that survive filtering. Without this, filtered
    % pixels appear "active" and get tau_active instead of tau_release.
    masked_input = polarity_map .* filter_mask;

    activity_indicator = single(masked_input ~= 0);
    activity_indicator(isnan(activity_indicator)) = 0;

    % Morphological dilation bridges sub-pixel gaps without the
    % amplitude amplification that Gaussian blur would cause.
    se = strel('disk', 1);
    activity_blurred = imdilate(activity_indicator, se);

    % ----------------------------------------------------------------
    % 3. Asymmetric attack-release envelope
    % ----------------------------------------------------------------
    tau_effective = tau_active .* activity_blurred + ...
                    surface_tau_release .* (1 - activity_blurred);

    % ----------------------------------------------------------------
    % 4. Coupled gain-decay coefficient
    % ----------------------------------------------------------------
    adaptive_gains = 1 - exp(-dt ./ tau_effective);

    % Inactive pixels: force alpha to zero (sample-and-hold)
    adaptive_gains(activity_blurred == 0) = 0;

    % Smooth the blending coefficient to soften spatial boundaries
    adaptive_gains = imgaussfilt(adaptive_gains, filter_sigma, ...
        'FilterSize', filter_size);

    % ----------------------------------------------------------------
    % 5. EMA surface update (first-order IIR)
    % ----------------------------------------------------------------
    time_surface_map_raw = adaptive_gains .* masked_input ...
        + (1 - adaptive_gains) .* time_surface_map_prev;

    % Suppress inactive pixels to eliminate trailing artifacts
    time_surface_map_raw = ...
        time_surface_map_raw .* single(activity_blurred > 0);

    % ----------------------------------------------------------------
    % 6. Divisive normalization (Carandini-Heeger)
    % ----------------------------------------------------------------
    magnitude = abs(time_surface_map_raw);

    activity_pool = imgaussfilt(magnitude, filter_sigma, ...
        'FilterSize', filter_size);
    counts_smooth = imgaussfilt(single(counts), filter_sigma, ...
        'FilterSize', filter_size);

    % Semi-saturation constant from median activity
    time_surface_map = zeros(size(masked_input));
    good_mask = (abs(masked_input) > 0);
    sigma = median(activity_pool(abs(activity_pool) > 0));

    time_surface_map(good_mask) = ...
        time_surface_map_raw(good_mask) ./ ...
        (sigma + counts_smooth(good_mask) .^ div_norm_exp);

    % ----------------------------------------------------------------
    % 7. Outlier rejection
    % ----------------------------------------------------------------
    mean_value_pos = mean(time_surface_map(time_surface_map > 0));
    mean_value_neg = mean(time_surface_map(time_surface_map < 0));
    std_value_pos  = std(time_surface_map(time_surface_map > 0));
    std_value_neg  = std(time_surface_map(time_surface_map < 0));

    pos_threshold = mean_value_pos + 4 * std_value_pos;
    neg_threshold = mean_value_neg - 4 * std_value_neg;

    % Clamp outliers to the median of their polarity
    med_pos = median(time_surface_map(time_surface_map > 0));
    med_neg = median(time_surface_map(time_surface_map < 0));

    time_surface_map(time_surface_map > pos_threshold) = med_pos;
    time_surface_map(time_surface_map < neg_threshold) = med_neg;

    time_surface_map_raw(time_surface_map > pos_threshold) = ...
        median(time_surface_map_raw(time_surface_map > 0));
    time_surface_map_raw(time_surface_map < neg_threshold) = ...
        median(time_surface_map_raw(time_surface_map < 0));

    % ----------------------------------------------------------------
    % 8. Tone mapping for display
    % ----------------------------------------------------------------
    % Bipolar surface mapped to [0, 1] via symmetric sigmoid.
    % Zero maps to 0.5 (mid-gray).
    normalized_output_frame = ...
        process.symmetricToneMappingNorm(time_surface_map, ...
        symmetric_tone_scale);

end