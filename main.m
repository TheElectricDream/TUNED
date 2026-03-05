clear;
clc;
close all;

%% Define processing range
% Define start and end time to process [seconds]
t_start_process = 80; 
t_end_process   = 100; 

%% Import events for inspection

% Set path to datasets
hdf5Path = ['/home/alexandercrain/Dropbox/Graduate Documents' ...
    '/Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5/'];
videoOutPath = '/home/alexandercrain/Videos/Research/';

% Set dataset name
fileName = 'recording_20260127_145247.hdf5';  % Jack W. (LED Cont)
%fileName = 'recording_20251029_131131.hdf5';  % EVOS - NOM - ROT
%fileName = 'recording_20251029_135047.hdf5';  % EVOS - SG - ROT
%fileName = 'recording_20251029_134602.hdf5';  % EVOS - DARK - ROT

% Load the data
tk = double(h5read([hdf5Path fileName], '/timestamp'));
xk = single(h5read([hdf5Path fileName], '/x'));
yk = single(h5read([hdf5Path fileName], '/y'));
pk = single(h5read([hdf5Path fileName], '/polarity'));

% Convert time to seconds
tk = (tk - tk(1))/1e6;

% Convert to single data type to use less memory
tk = single(tk);

% Find indices within the valid range
valid_idx = tk >= t_start_process & tk <= t_end_process;

% Filter the data vectors
tk = tk(valid_idx);
xk = xk(valid_idx);
yk = yk(valid_idx);
pk = pk(valid_idx);

% Shift time to start at 0 for the new window
% This ensures your frame loop starts correctly at frame 1
tk = tk - t_start_process; 

% Clear unused variables for memory
clearvars valid_idx;

%% Initialize all tunable parameters for the algorithms
% GENERAL PARAMETERS
% ------------------
% Set the image size
imgSz                       = [640, 480]; 

% Set the time interval to accumulate over
t_interval                  = 0.33;  % [s]
t_total                     = max(tk);  % [s]
frame_total                 = floor(t_total/t_interval);

% Boolean controls
filter_output_image         = false;

% INTER-EVENT-INTERVAL
% --------------------
% Persistent inter-event-interval map (EMA) parameter
iei_alpha                   = 0.8;     

% Initialize EMA-IEI 
iei_map                     = zeros(imgSz);

% COHERENCE PARAMETERS
% --------------------
% Define coherence parameters - Must be tuned based on t_interval
coh_params.r_s                         = 30/imgSz(1);  % spatial radius [pixels norm]
coh_params.trace_threshold             = 1.3;
coh_params.persistence_threshold       = 0.0002;
coh_params.coherence_threshold         = 0.06;
coh_params.similarity_threshold        = 0.5;

% Initialize mask for filter
filter_mask                 = ones(imgSz);

% Boolean controls
filter_by_coherence         = true;

% ADAPTIVE LOCAL TIME-SURFACE PARAMETERS
% --------------------------------------
% Set adaptive local time-surface parameters
alts_params.dt                   = t_interval;
alts_params.filter_size          = 11;
alts_params.filter_sigma         = 9.0;
alts_params.surface_tau_release  = 3.0;
alts_params.div_norm_exp         = 1.0;
alts_params.symmetric_tone_scale = 3.0;
alts_activity_score.mean         = zeros(frame_total, 1);
alts_activity_score.median       = zeros(frame_total, 1);
alts_activity_score.std          = zeros(frame_total, 1);

% Create a cell array to store per-frame data (preallocate for frame_total)
alts_frame_storage      = cell(frame_total,1);

% TIME SURFACE PARAMETERS
% -----------------------
% Initialize map
ts_t_map = -inf(imgSz); 

% Decay constant for the visual 
ts_time_constant = 0.05;  % [seconds]

% SPEED INVARIENT TIME SURFACE PARAMETERS
% ---------------------------------------
% REF: https://arxiv.org/pdf/1903.11332
% ---------------------------------------
% Initialize map
sits_t_map = zeros(imgSz);

% Radius of neighbourhood
sits_R = 3;

% ADAPTIVE GLOBAL DECAY TIME SURFACE PARAMETERS
% ---------------------------------------------
% REF: https://ieeexplore.ieee.org/document/10205486/
% ---------------------------------------------
% Initialize map
agd_surface = zeros(imgSz); 

% State structure to hold history
agd_state.last_t_map = zeros(imgSz);  % Stores timestamp of last event
agd_state.activity = 0;  % Initial Activity level
agd_state.last_update_time = t_start_process; 

% Tuning parameters
% Note: for unfiltered events with large spikes, an aggressive smoothing 
% factor is needed. Otherwise, the scene rests anytime it sees a hot pixel
agd_params.alpha   = 0.001;   % Smoothing factor for activity 
agd_params.K       = 5000000.0;  % Scaling factor (Controls "memory length")
agd_activity_store = zeros(length(frame_total),1);

% MOTION-ENCODED TIME-SURFACE (METS) PARAMETERS
% ----------------------------------------------
% REF: Xu et al., "METS: Motion-Encoded Time-Surface for Event-Based
%      High-Speed Pose Tracking," IJCV, vol. 133, pp. 4401-4419, 2025.
%      DOI: 10.1007/s11263-025-02379-6
% ----------------------------------------------
% State: polarity-separated timestamp maps (Eq. 4)
mets_state.t_last_pos = zeros(imgSz);   % Last event timestamp, positive polarity
mets_state.t_last_neg = zeros(imgSz);   % Last event timestamp, negative polarity
mets_state.t_last_any = zeros(imgSz);   % Last event timestamp, either polarity
mets_state.p_last     = zeros(imgSz);   % Polarity of last event at each pixel

% Parameters (Section 3.3 of the paper — identical to their defaults)
mets_params.R    = 4;    % Observation window half-size [pixels]
mets_params.n    = 3;    % Velocity estimation range [pixels]
mets_params.d    = 5;    % Decay step [pixels]
mets_params.d_th = 8;    % Decay distance threshold [pixels]

%% Initialize all tunable parameters for feature detection algorithms
% ARC* CORNER DETECTION PARAMETERS
% ---------------------------------
% Reference: Alzugaray & Chli, IEEE RA-L 2018
%
% kappa:  S* filter refractory period [seconds]. Events at the same pixel
%         arriving within kappa of each other are considered redundant and
%         skipped. Higher values = more aggressive filtering. The paper
%         uses 50 ms. For the EVOS dataset with slow spacecraft motion,
%         values of 10-50 ms are reasonable.
%
% radii:  Bresenham circle radii for the dual-ring test. Both rings must
%         independently classify the event as a corner.
%
% arc_bounds: {Mx1} cell of [Lmin, Lmax] per radius. These control the
%             angular range of detectable corners. Lmin prevents noise
%             from triggering tiny arcs. Lmax limits the maximum arc
%             angle (preventing edges from being classified as corners).
%             With 16 elements (r=3): [3,6] → approx 67°–135° arc range.

arc_star_params.kappa      = 0.050;         % S* refractory period [s]
arc_star_params.radii      = [3, 4];        % Circle radii
arc_star_params.arc_bounds = {[4, 6], [5, 8]};

% Boolean controls
detect_arc_star_corners    = false;   % Enable/disable Arc* detection
overlay_arc_star_corners   = false;   % Enable/disable corner overlay on video

% Marker appearance
marker_opts.radius    = 2;               % Cross arm length [pixels]
marker_opts.color_pos = uint8([0 255 0]);    % Positive polarity (green)
marker_opts.color_neg = uint8([255 0 255]);  % Negative polarity (magenta)
marker_opts.thickness = 1;               % Cross line width [pixels]

% Storage for detected corner events across all frames
corner_events_storage = cell(frame_total, 1);

% HARRIS-CV CORNER DETECTION PARAMETERS
% --------------------------------------
% Boolean: enable/disable Harris-CV detection
detect_harris_cv_corners       = true;

% Boolean: burn markers onto the ATS output frame
overlay_harris_cv_corners      = true;

% Harris detector parameters (see +features/detectHarrisCV.m)
harris_cv_params.k             = 0.1;     % Harris sensitivity [0.04 - 0.15]
harris_cv_params.sigma_smooth  = 1.5;      % Gaussian integration window [px]
harris_cv_params.threshold     = 0.0001;   % Fraction of max(R) to keep
harris_cv_params.nms_radius    = 8;        % Non-max suppression radius [px]
harris_cv_params.max_corners   = 150;      % Max corners per frame
harris_cv_params.border        = 10;       % Border exclusion [px]

% Storage for per-frame corners (optional, for post-processing)
harris_cv_storage = cell(frame_total, 1);

% CV DESCRIPTOR PARAMETERS
harris_cv_desc_params.half_size      = 5;      % 11×11 = 121-dim descriptor
harris_cv_desc_params.sigma_weight   = 3.0;    % Gaussian centre-weighting
harris_cv_desc_params.min_observed   = 0.25;   % Min observed fraction

% MATCHING PARAMETERS
harris_cv_match_params.ncc_threshold = 0.5;    % Min NCC to accept match
harris_cv_match_params.ratio_test    = 1.0;    % Lowe's ratio (IJCV 2004)
harris_cv_match_params.max_distance  = 5;     % Max displacement [px]
harris_cv_match_params.mutual        = true;   % Symmetric matching

% Display control
overlay_harris_cv_matches            = true;

% Previous-frame state (persists across iterations)
desc_prev_frame                      = [];
corners_prev_frame                   = [];
corners_prev_frame_for_display       = [];

% Per-frame match storage
harris_cv_match_storage = cell(frame_total, 1);

% DEBUG VIEWER
% Enable this to step through frames interactively.
% Set to false for batch processing / video recording.
debug_harris_cv_matches         = true;

% Storage for the previous frame's ATS output (sensor orientation,
% BEFORE transpose). This is needed so the viewer can show the
% actual image content of frame N-1 alongside frame N.
debug_prev_frame                = [];

% LEVEL 3: MULTI-CHANNEL IEI DESCRIPTOR PARAMETERS
% --------------------------------------------------
% Patch geometry (same for all channels)
harris_cv_desc_params.half_size      = 5;      % 11×11 per channel
harris_cv_desc_params.sigma_weight   = 3.0;    % Gaussian centre-weight
harris_cv_desc_params.min_observed   = 0.25;   % Min obs fraction

% Channels to include in the descriptor.
% Each must correspond to a field in the maps struct built in the loop.
%
%   'cv'       — CV map (σ/μ of local IEI). Regularity pattern.
%   'iei_mean' — Persistent EMA IEI map. Absolute timing / speed.
%   'iei_std'  — Current-frame IEI std. Temporal variability.
%   'harris_r' — Harris response on CV map. Corner shape in CV space.
%   'density'  — Trace map (spatial density). Activity level.
%
harris_cv_desc_params.channels = {'cv', 'iei_mean', 'iei_std', 'harris_r', 'density'};

% Per-channel weights (optional). Channels with higher weight
% contribute more to the final descriptor after normalisation.
% All default to 1.0 if not specified.
%
% Tuning strategy:
%   - Start with all 1.0 (equal contribution)
%   - If matches are still ambiguous along edges, INCREASE iei_mean
%     weight (breaks speed degeneracy)
%   - If corners at junctions match well but edges don't, INCREASE
%     harris_r weight (encodes corner geometry)
%   - If you get false matches in low-activity regions, ADD density
harris_cv_desc_params.channel_weights = struct( ...
    'cv',       1.0, ...
    'iei_mean', 1.0, ...
    'iei_std',  1.0, ...
    'harris_r', 1.0, ...
    'density', 1.0);

%% Initialize all figure code for video output
% Indicate which videos should be saved
cohOut = false;
atsOut = true;

% Initialize the videos
[hFigs, hAxs, hImgs, videoWriters] = plot.initializeEventVideos(cohOut,...
    atsOut, imgSz, videoOutPath);

%% Initialize data storage and perform data optimizations
% Identify number of events
current_idx     = 1;
n_events        = length(tk);
frame_time      = zeros(frame_total, 1);

% Initialize per-pixel timestamp tracking
last_event_timestamp    = zeros(imgSz);
norm_trace_map_prev     = zeros(imgSz);
time_surface_map_prev   = zeros(imgSz);

% Initialize polarity-separated Surface of Active Events (S*)
% for the Arc* corner detector. Persists across frames.
sae_state.t_last    = zeros(imgSz);  % Last event time (any polarity)
sae_state.p_last    = zeros(imgSz);  % Last event polarity at each pixel
sae_state.t_ref_pos = zeros(imgSz);  % Reference time, positive polarity
sae_state.t_ref_neg = zeros(imgSz);  % Reference time, negative polarity

%% Data processing starting point
% Loop through the figures to capture each frame
for frameIndex = 1:frame_total  

    % Start loop timer
    tic;
    
    % Increment the interval
    t_range_c = (frameIndex - 1) * t_interval;
    t_range_n = (t_range_c+t_interval);

    % Slice the events to a valid range
    [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    process.sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, current_idx);

    % Confirm the presence of valid events in the packet
    % If no events are present, we skip this frame
    if isempty(t_valid)
        
        fprintf('There are no events in this slice, skipping... \n');
        continue;
        
    end   

    % ---------------------- EVENT PREPERATION------------------------%
    % ----------------------------------------------------------------%
    
    % Convert 2D subscripts (x,y) to 1D linear indices
    % Imagine you are a post-man with a disorganized stack of letters. 
    % Instead of dealing with letters for 3rd Avenue, 5th street, the 
    % sub2ind function assigns and "ID" for each house. So going forward, 
    % (3,5) might just be "House #1".Programatically this just means that
    % (1,1) is "1", (1,2) is "2". So you will have a list at the end which
    % is of size x*y. In this case with a frame of size 640 by 480, the 
    % MAXIMUM size of the list is 307,200. However the practical size of 
    % the list will change as not every pixel is active.
    linear_idx = sub2ind(imgSz, x_valid, y_valid);
    
    % So now that we have a linear index, we sort them by "ID". So
    % ultimately what we get here is a sorted list of event times where we
    % have all event time differences for each pixel coordinate groupped. 
    % So the time intervals may not actually be monotonic. So
    % sorted_idx is basically the list of "houses" grouped together.
    [sorted_idx, sort_order] = sort(linear_idx);
    sorted_t = t_valid(sort_order);
    sorted_x = x_valid(sort_order);
    sorted_y = y_valid(sort_order);
    sorted_p = p_valid(sort_order);

    % Ensure polarity is -1 and 1 (if it's 0 and 1)
    p_signed = double(sorted_p);
    p_signed(p_signed == 0) = -1;

    % Reset the frames for the current loop
    counts = zeros(imgSz);
    
    % With the sorted index list, aka the houses, we can extract where
    % in that list each "house" starts. That would be the "pos" output.
    % The unique_idx list is the complete list of all houses. 
    [unique_idx, pos, ~] = unique(sorted_idx);
           
    % Now we need to find where each "house" starts and ends. The start
    % is easy because we get it from "pos". The end can be inferred
    % from the start because the start of the NEXT "house" group is
    % always one more than the end of the previous group. So,
    % pos(2:end)-1 gives the list of ends. We do not know when the last
    % group will end, but we do know that it MUST end by the end of the
    % dataset. So this is the length(sorted_idx).
    group_ends = [pos(2:end)-1; length(sorted_idx)];

    % This is a bit of cheeky MATLAB code. Because "unique_idx" is a
    % linear index list, MATLAB knows to map it to the size of the 2D
    % array "counts". So although group_end and pos are vectors, the
    % output is a 640x480 array. Additionally, the ending position of
    % each "house" minus the starting position (i.e. the difference)
    % gives you the total number of "things" assigned to that "ID". 
    counts(unique_idx) = group_ends - pos + 1;

    % ------------------------- STATISTICS -------------------------------%
    % --------------------------------------------------------------------%

    [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
        stats.computeNeighborhoodStats(sorted_t, unique_idx, pos, ...
        group_ends, imgSz);

    % Update a persistant map of the inter-event interval so that sparse
    % data is retained
    new_obs_mask = (t_mean_diff > 0); 
    iei_map(new_obs_mask) = (1 - iei_alpha) .* iei_map(new_obs_mask) +...
        iei_alpha .* t_mean_diff(new_obs_mask);

    % ---------------------- EVENT COHERENCE -----------------------------%
    % --------------------------------------------------------------------%

    % Choose whether or not to generate a filter
    if filter_by_coherence == true

        [norm_trace_map, norm_similarity_map, norm_persist_map,...
            filtered_coherence_map] = coherence.computeCoherenceMask(sorted_x,...
            sorted_y, sorted_t, imgSz, t_interval, unique_idx, pos, ...
            group_ends, coh_params, frameIndex, norm_trace_map_prev, t_mean_diff);

        % Set any retention variables
        norm_trace_map_prev = norm_trace_map;

        % Create the filter
        filter_mask = filtered_coherence_map;
        filter_mask(isnan(filter_mask)) = 0;
        filter_mask = single(imgaussfilt(single(filter_mask), 5.0, "FilterSize", 9));
        filter_mask(filter_mask<coh_params.coherence_threshold) = 0;
    else

        % Use a unity mask instead
        filter_mask = ones(imgSz);

    end

    % ---------------------- ARC* CORNER DETECTION ----------------------%
    % -------------------------------------------------------------------%
    
    % if detect_arc_star_corners
    %     % Arc* processes events in temporal order (not pixel-sorted order)
    %     % so we pass the original slice outputs directly.
    %     [corner_events_frame, sae_state] = features.detectArcStarCorners(...
    %         x_valid, y_valid, t_valid, p_valid, imgSz, ...
    %         sae_state, arc_star_params);
    % 
    %     % Store corner events for this frame
    %     corner_events_storage{frameIndex} = corner_events_frame;
    % 
    %     % Report detection rate
    %     if ~isempty(corner_events_frame)
    %         n_events_frame = numel(t_valid);
    %         n_corners_frame = size(corner_events_frame, 1);
    %         corner_pct = 100 * n_corners_frame / n_events_frame;
    %         fprintf('  Arc*: %d corners / %d events (%.1f%%)\n', ...
    %             n_corners_frame, n_events_frame, corner_pct);
    %     end
    % end

    % ------------------- HARRIS-CV DETECTION + MATCHING ----------------%
    % -------------------------------------------------------------------%

    if detect_harris_cv_corners

        % --- Step 1: CV map from persistent IEI (UNCHANGED) ---
        [~, cv_map_for_harris, ~] = coherence.findSimilarities( ...
            sorted_x, sorted_y, t_mean_diff.*filter_mask, imgSz, 8);
        cv_map_for_harris(isnan(cv_map_for_harris)) = 0;
        cv_map_for_harris(isinf(cv_map_for_harris)) = 0;

        % --- Step 2: Detect corners (UNCHANGED) ---
        [harris_corners, harris_R, harris_corners_sub] = ...
            features.detectHarrisCV(cv_map_for_harris, ...
            'k',            harris_cv_params.k, ...
            'sigma_smooth', harris_cv_params.sigma_smooth, ...
            'threshold',    harris_cv_params.threshold, ...
            'nms_radius',   harris_cv_params.nms_radius, ...
            'max_corners',  harris_cv_params.max_corners, ...
            'border',       harris_cv_params.border);

        harris_cv_storage{frameIndex} = harris_corners;

        % --- Step 3: Extract MULTI-CHANNEL descriptors (CHANGED) ---
        %
        %  Build the maps struct from data already computed earlier
        %  in this loop iteration. All of these variables exist by
        %  this point in the pipeline:
        %
        %    cv_map_for_harris  — computed in Step 1 above
        %    iei_map            — persistent EMA, updated before coherence
        %    t_std_diff         — from stats.computeNeighborhoodStats
        %    harris_R           — from Step 2 above
        %    norm_trace_map     — from coherence.computeCoherenceMask
        %
        if ~isempty(harris_corners_sub)

            % Assemble channel maps
            iei_desc_maps = struct();
            iei_desc_maps.cv       = cv_map_for_harris;
            iei_desc_maps.iei_mean = t_mean_diff.*filter_mask;
            iei_desc_maps.harris_r = harris_R;

            % For iei_std: use t_std_diff (current-frame IEI variability)
            % from computeNeighborhoodStats. If you named it differently,
            % adjust here. Replace NaN with 0 for consistency.
            iei_std_map = t_std_diff;
            iei_std_map(isnan(iei_std_map)) = 0;
            iei_desc_maps.iei_std = iei_std_map;

            % Optional: add density channel if enabled
            if ismember('density', harris_cv_desc_params.channels)
                iei_desc_maps.density = norm_trace_map;
            end

            % Extract multi-channel descriptors
            [desc_curr, corners_curr_valid] = ...
                features.extractIEIDescriptors(iei_desc_maps, ...
                harris_corners_sub, ...
                'half_size',       harris_cv_desc_params.half_size, ...
                'sigma_weight',    harris_cv_desc_params.sigma_weight, ...
                'min_observed',    harris_cv_desc_params.min_observed, ...
                'channels',        harris_cv_desc_params.channels, ...
                'channel_weights', harris_cv_desc_params.channel_weights);
        else
            desc_curr = [];
            corners_curr_valid = zeros(0, 2);
        end

        % --- Step 4: Match (UNCHANGED — works on any dimension) ---
        harris_cv_matches = zeros(0, 2);
        harris_cv_scores  = zeros(0, 1);

        if ~isempty(desc_curr) && ~isempty(desc_prev_frame)
            [harris_cv_matches, harris_cv_scores] = ...
                features.matchCVDescriptors( ...
                desc_curr, corners_curr_valid, ...
                desc_prev_frame, corners_prev_frame, ...
                'ncc_threshold', harris_cv_match_params.ncc_threshold, ...
                'ratio_test',    harris_cv_match_params.ratio_test, ...
                'max_distance',  harris_cv_match_params.max_distance, ...
                'mutual',        harris_cv_match_params.mutual);
        end

        % Store for post-processing
        harris_cv_match_storage{frameIndex} = struct( ...
            'matches', harris_cv_matches, ...
            'scores',  harris_cv_scores, ...
            'corners_curr', corners_curr_valid, ...
            'corners_prev', corners_prev_frame);

        % Save for display BEFORE overwriting
        corners_prev_frame_for_display = corners_prev_frame;

        % Update previous-frame state
        desc_prev_frame    = desc_curr;
        corners_prev_frame = corners_curr_valid;

        % Report (now shows descriptor dimension)
        n_det = size(corners_curr_valid, 1);
        n_match = size(harris_cv_matches, 1);
        d_dim = size(desc_curr, 2);
        if n_match > 0
            fprintf('  Harris-CV L3 [%d-dim]: %d det, %d matched (NCC: %.3f)\n', ...
                d_dim, n_det, n_match, mean(harris_cv_scores));
        elseif n_det > 0
            fprintf('  Harris-CV L3 [%d-dim]: %d det, 0 matched\n', ...
                d_dim, n_det);
        end

    end

    % -------------- ADAPTIVE LOCAL TIME-SURFACE UPDATE ------------------%
    % --------------------------------------------------------------------%
    
    % Accumulate polarity into a 2D grid
    % If multiple events land on one pixel, we sum their polarities (e.g., +1 +1 -1 = +1)
    polarity_map = accumarray([sorted_x, sorted_y], p_signed, imgSz, @sum, 0);

    [normalized_output_frame, time_surface_map_raw, tau_filtered, adaptive_gains] = ...
    accumulator.localAdaptiveTimeSurface(iei_map,...
    time_surface_map_prev, alts_params, filter_mask, polarity_map, counts);

    % Set any retention variables
    time_surface_map_prev = time_surface_map_raw;

    % -------------------------- DATA EXPORT -----------------------------%
    % --------------------------------------------------------------------%
    
    % Store the processed frames for later use
    alts_frame_storage{frameIndex} = normalized_output_frame;

    % Store the adaptive map score
    alts_activity_score.mean(frameIndex) = mean(adaptive_gains(abs(adaptive_gains)>0));
    alts_activity_score.median(frameIndex) = median(adaptive_gains(abs(adaptive_gains)>0));
    alts_activity_score.std(frameIndex) = std(adaptive_gains(abs(adaptive_gains)>0));

    % ------------------------ EXPORTING VIDEO ---------------------------%
    % --------------------------------------------------------------------%

    if debug_harris_cv_matches

        plot.debugMatchViewer( ...
            debug_prev_frame, ...               % previous ATS frame
            normalized_output_frame, ...        % current ATS frame
            corners_prev_frame_for_display, ...             % prev corners
            corners_curr_valid, ...             % curr corners
            harris_cv_matches, ...              % match pairs
            harris_cv_scores, ...               % NCC scores
            frameIndex);

        % Store current frame for next iteration's "previous"
        debug_prev_frame = normalized_output_frame;

    end    
    
    % Store the processed frames for later use
    alts_frame_storage{frameIndex} = normalized_output_frame;

    % Convert to uint8 (0-255 range)
    grayscale_normalized_output_frame = ...
        uint8(normalized_output_frame .* 255);

    % Apply a Gaussian filter to help smooth out the final image
    if filter_output_image == true 
         grayscale_normalized_output_frame = ...
             imgaussfilt(grayscale_normalized_output_frame, ...
             3.0,"FilterSize",3); 
    end

    % Log processing time
    frame_time(frameIndex) = toc;

    % Print progress
    stats.printPercentComplete(frameIndex, frame_total, frame_time(frameIndex));

end

% Close the video writer
for videosIdx = 1:length(videoWriters)
    close(videoWriters{videosIdx});
end


