clear;
clc;
close all;

%% Processing pipeline reminder
% To import AEDAT4 -->
% /home/alexandercrain/Repositories/CNN/import/importAEDAT4toHDF5.py

% All datasets --> /home/alexandercrain/Dropbox/Graduate Documents/
% Doctor of Philosophy/Thesis Research/Datasets/SPOT

% Datasets for MATLAB --> /home/alexandercrain/Dropbox/Graduate Documents/
% Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5

% Notes:
% - Look into kriging interpolation

%% Define processing range
% Define start and end time to process [seconds]
t_start_process = 80; 
t_end_process   = 1000; 

%% Import events for inspection

% Set path to datasets
hdf5_path = ['/home/alexandercrain/Dropbox/Graduate Documents' ...
    '/Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5/'];

% Set dataset name
%file_name = 'recording_20260127_145247.hdf5';  % Jack W. (LED Cont)
file_name = 'recording_20251029_131131.hdf5';  % EVOS - NOM - ROT
%file_name = 'recording_20251029_135047.hdf5';  % EVOS - SG - ROT
%file_name = 'recording_20251029_134602.hdf5';  % EVOS - DARK - ROT

% Load the data
tk = double(h5read([hdf5_path file_name], '/timestamp'));
xk = single(h5read([hdf5_path file_name], '/x'));
yk = single(h5read([hdf5_path file_name], '/y'));
pk = single(h5read([hdf5_path file_name], '/polarity'));

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

%% Initialize all tunable parameters for frame generation & filtering
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
coh_params.coherence_threshold         = 0.0001;
coh_params.similarity_threshold        = 0.7;

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
harris_cv_params.k             = 0.04;     % Harris sensitivity [0.04 - 0.15]
harris_cv_params.sigma_smooth  = 1.5;      % Gaussian integration window [px]
harris_cv_params.threshold     = 0.0001;     % Fraction of max(R) to keep
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
harris_cv_match_params.max_distance  = 10;     % Max displacement [px]
harris_cv_match_params.mutual        = false;   % Symmetric matching

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

%% Initialize all figure code for video output
% Indicate which videos should be saved
cohOut = false;
atsOut = true;

% Initialize the videos
[hFigs, hAxs, hImgs, videoWriters] = plot.initializeEventVideos(cohOut,...
    atsOut, imgSz);

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

    % ------------------- HARRIS-CV CORNER DETECTION --------------------%
    % -------------------------------------------------------------------%

    % if detect_harris_cv_corners
    % 
    %     % Compute the CV map from the persistent IEI map.
    %     % We call findSimilarities directly rather than extracting
    %     % the CV map from computeCoherenceMask, because:
    %     %   (a) computeCoherenceMask uses t_mean_diff (current-frame IEI)
    %     %       for coherence scoring, but we want the smoothed iei_map
    %     %       for corner detection (more stable statistics).
    %     %   (b) This avoids modifying the coherence pipeline.
    %     %
    %     % The radius here controls the spatial scale of the CV
    %     % computation. Using a slightly larger radius (6-10) than the
    %     % default (4) gives smoother CV fields that produce cleaner
    %     % Harris gradients.
    %     [~, cv_map_for_harris, ~] = coherence.findSimilarities( ...
    %         sorted_x, sorted_y, t_mean_diff.*filter_mask, imgSz, 8);
    % 
    %     % Replace NaN (unobserved pixels) with 0
    %     cv_map_for_harris(isnan(cv_map_for_harris)) = 0;
    %     cv_map_for_harris(isinf(cv_map_for_harris)) = 0;
    % 
    %     % Run Harris on the CV map
    %     [harris_corners, harris_R, harris_corners_sub] = ...
    %         features.detectHarrisCV(cv_map_for_harris, ...
    %         'k',            harris_cv_params.k, ...
    %         'sigma_smooth', harris_cv_params.sigma_smooth, ...
    %         'threshold',    harris_cv_params.threshold, ...
    %         'nms_radius',   harris_cv_params.nms_radius, ...
    %         'max_corners',  harris_cv_params.max_corners, ...
    %         'border',       harris_cv_params.border);
    % 
    %     % Store for post-processing (optional)
    %     harris_cv_storage{frameIndex} = harris_corners;
    % 
    %     % Report
    %     if ~isempty(harris_corners)
    %         fprintf('  Harris-CV: %d corners detected\n', ...
    %             size(harris_corners, 1));
    %     end
    % 
    % end

    % ------------------- HARRIS-CV DETECTION + MATCHING ----------------%
    % -------------------------------------------------------------------%

    if detect_harris_cv_corners

        % --- CV map from persistent IEI ---
        [~, cv_map_for_harris, ~] = coherence.findSimilarities( ...
            sorted_x, sorted_y, t_mean_diff.*filter_mask, imgSz, 8);
        cv_map_for_harris(isnan(cv_map_for_harris)) = 0;
        cv_map_for_harris(isinf(cv_map_for_harris)) = 0;

        % --- Detect corners (Level 1) ---
        [harris_corners, harris_R, harris_corners_sub] = ...
            features.detectHarrisCV(cv_map_for_harris, ...
            'k',            harris_cv_params.k, ...
            'sigma_smooth', harris_cv_params.sigma_smooth, ...
            'threshold',    harris_cv_params.threshold, ...
            'nms_radius',   harris_cv_params.nms_radius, ...
            'max_corners',  harris_cv_params.max_corners, ...
            'border',       harris_cv_params.border);

        harris_cv_storage{frameIndex} = harris_corners;

        % --- Extract descriptors (Level 2) ---
        if ~isempty(harris_corners_sub)
            [desc_curr, corners_curr_valid] = ...
                features.extractCVDescriptors(cv_map_for_harris, ...
                harris_corners_sub, ...
                'half_size',    harris_cv_desc_params.half_size, ...
                'sigma_weight', harris_cv_desc_params.sigma_weight, ...
                'min_observed', harris_cv_desc_params.min_observed);
        else
            desc_curr = [];
            corners_curr_valid = zeros(0, 2);
        end

        % --- Match to previous frame (Level 2) ---
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

        % Save previous-frame state for display BEFORE overwriting
        corners_prev_frame_for_display = corners_prev_frame;

        % Update previous-frame state
        desc_prev_frame    = desc_curr;
        corners_prev_frame = corners_curr_valid;

        % Report
        n_det = size(corners_curr_valid, 1);
        n_match = size(harris_cv_matches, 1);
        if n_match > 0
            fprintf('  Harris-CV: %d det, %d matched (NCC: %.3f)\n', ...
                n_det, n_match, mean(harris_cv_scores));
        elseif n_det > 0
            fprintf('  Harris-CV: %d det, 0 matched\n', n_det);
        end

        if ~isempty(harris_corners)
            fprintf('  Harris-CV: %d corners detected\n', ...
                size(harris_corners, 1));
        end

    end

    % -------------- ADAPTIVE LOCAL TIME-SURFACE UPDATE ------------------%
    % --------------------------------------------------------------------%
    
    % Accumulate polarity into a 2D grid
    % If multiple events land on one pixel, we sum their polarities
    % (e.g., +1 +1 -1 = +1)
    polarity_map = accumarray([sorted_x, sorted_y], p_signed, imgSz, @sum, 0);

    [normalized_output_frame, time_surface_map_raw, tau_filtered, adaptive_gains] = ...
    accumulator.localAdaptiveTimeSurface(iei_map,...
    time_surface_map_prev, alts_params, filter_mask, polarity_map, counts);

    % Set any retention variables
    time_surface_map_prev = time_surface_map_raw;

    % Store the adaptive map score
    alts_activity_score.mean(frameIndex) = ...
        mean(adaptive_gains(abs(adaptive_gains)>0));
    alts_activity_score.median(frameIndex) = ...
        median(adaptive_gains(abs(adaptive_gains)>0));
    alts_activity_score.std(frameIndex) = ...
        std(adaptive_gains(abs(adaptive_gains)>0));

    % ----------------- NUNES GLOBAL ADAPTIVE ACCUMULATION----------------%
    % --------------------------------------------------------------------%
    
    % % Run the AGD algorithm
    % [agd_surface, agd_state, ~] = accumulator.adaptiveGlobalDecay(agd_surface,...
    %     sorted_x, sorted_y, sorted_t, agd_state, agd_params);
    % 
    % % Store the surface into the standard normalized frame
    % normalized_output_frame = agd_surface;
    % 
    % % Store activity data for later inspection
    % agd_activity_store(frameIndex) = agd_state.activity;

    % --------------------- TIME-SURFACE ACCUMULATION --------------------%
    % --------------------------------------------------------------------%
    
    % % Run the normal time-surface accumulation algorithm
    % [ts_t_map, normalized_output_frame] = accumulator.timeSurface(ts_t_map,...
    %  sorted_x, sorted_y, sorted_t, imgSz, ts_time_constant);

    % ---------- SPEED INVARIENT TIME-SURFACE ACCUMULATION ---------------%
    % --------------------------------------------------------------------%
    
    % % Run the speed invarient time-surface accumulation algorithm
    % [sits_t_map, normalized_output_frame] = ...
    %     accumulator.speedInvariantTimeSurface(sits_t_map, sorted_x,...
    %     sorted_y, sits_R);

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

    % if atsOut
    %     % Transpose to display orientation (H x W)
    %     display_frame = grayscale_normalized_output_frame';
    % 
    %     % Overlay Arc* corner markers if enabled
    %     if overlay_arc_star_corners && detect_arc_star_corners ...
    %             && ~isempty(corner_events_storage{frameIndex})
    %         if exist('marker_opts', 'var')
    %             display_frame = plot.overlayCornerMarkers(...
    %                 display_frame, corner_events_storage{frameIndex}, ...
    %                 marker_opts);
    %         else
    %             display_frame = plot.overlayCornerMarkers(...
    %                 display_frame, corner_events_storage{frameIndex});
    %         end
    %     end
    % 
    %     % Overlay Harris-CV corner markers if enabled
    %     if overlay_harris_cv_corners && detect_harris_cv_corners ...
    %             && exist('harris_corners_sub', 'var') ...
    %             && ~isempty(harris_corners_sub)
    %         display_frame = plot.overlayHarrisCVMarkers(...
    %             display_frame, harris_corners_sub);
    %     end
    % 
    %     % Overlay Harris-CV match lines
    %     if overlay_harris_cv_matches && detect_harris_cv_corners ...
    %             && exist('harris_cv_matches', 'var') ...
    %             && ~isempty(harris_cv_matches) ...
    %             && ~isempty(corners_prev_frame_for_display)
    %         display_frame = plot.overlayMatchedCorners( ...
    %             display_frame, ...
    %             corners_curr_valid, ...
    %             corners_prev_frame_for_display, ...
    %             harris_cv_matches, ...
    %             harris_cv_scores);
    %     end
    % 
    %     % Capture the frame for the video writer
    %     set(hImgs{1}, 'CData', display_frame);
    %     set(hImgs{1}, 'AlphaData', ...
    %         ~isnan(grayscale_normalized_output_frame'));
    %     set(hAxs{1}, 'Visible','off');
    %     colormap(gray);
    %     clim([0 255]);
    %     writeVideo(videoWriters{1}, display_frame);
    % 
    %     frameOutputFolder = '/home/alexandercrain/Videos/output_frames';
    %     if ~exist(frameOutputFolder, 'dir')
    %         mkdir(frameOutputFolder);
    %     end
    %     fname = fullfile(frameOutputFolder, sprintf('frame_%05d.png', frameIndex));
    %     imwrite(display_frame, fname);
    % end

    % Print progress
    stats.printPercentComplete(frameIndex, frame_total, frame_time(frameIndex));

end

% Concatenate all detected corner events into a single array
if detect_arc_star_corners
    all_corner_events = vertcat(corner_events_storage{:});
    fprintf('\nArc* total: %d corner events detected.\n', ...
        size(all_corner_events, 1));
    % Columns: [x, y, t, polarity]
end

if detect_harris_cv_corners
    total_matches = 0;
    all_ncc = [];
    for fi = 1:frame_total
        if ~isempty(harris_cv_match_storage{fi})
            ms = harris_cv_match_storage{fi};
            total_matches = total_matches + size(ms.matches, 1);
            all_ncc = [all_ncc; ms.scores]; %#ok<AGROW>
        end
    end
    fprintf('\n=== Harris-CV Level 2 Summary ===\n');
    fprintf('Total matches: %d across %d frames (%.1f/frame)\n', ...
        total_matches, frame_total, total_matches / frame_total);
    if ~isempty(all_ncc)
        fprintf('NCC: mean=%.3f  std=%.3f  [%.3f, %.3f]\n', ...
            mean(all_ncc), std(all_ncc), min(all_ncc), max(all_ncc));
    end
end



% Close the video writer
for videosIdx = 1:length(videoWriters)
    close(videoWriters{videosIdx});
end


