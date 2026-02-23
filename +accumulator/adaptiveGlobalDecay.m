function [S, state, normalized_output_frame] = adaptiveGlobalDecay(S, x_list, y_list, t_list, state, params)
% ADAPTIVEGLOBALDECAY Update surface using Global Event Activity.
%
%   [S, state] = adaptiveGlobalDecay(S, x, y, t, p, state, params)
%
%   Inputs:
%     S      - Current surface (H x W matrix)
%     x,y,t  - Event vectors (spatial coords, timestamps)
%     state  - Struct containing persistent variables (last_t, activity_state)
%     params - Struct with tuning parameters (alpha, beta, K, etc.)
%
%   Ref: Nunes et al., "Adaptive Global Decay Process for Event Cameras", CVPR 2023

    if isempty(t_list)
        return;
    end

    % 1. Determine the time delta for this batch
    t_batch_end = t_list(end);
    dt_batch = t_batch_end - state.last_update_time;

    % Handle duplicates and negative time
    % If dt is zero (duplicate packet) or negative (out of order), 
    % force it to be at least 1 microsecond (1e-6).
    if dt_batch < 1e-6
        fprintf('\ndt_batch is zero, setting to 1e-6\n')
    end
    dt_batch = max(dt_batch, 1e-6);

    % 2. Calculate Global Activity A(t)
    % N_events / (Total Pixels * Time) or just N_events / Time
    % The paper often normalizes by resolution.
    [H, W] = size(S);
    current_rate = numel(t_list) / dt_batch; 
    
    % Normalize rate (optional, helps keep params K stable across resolutions)
    % current_rate = current_rate / (H * W); 

    % Smooth the activity (Leaky Integrator)
    state.activity = (1 - params.alpha) * state.activity + params.alpha * current_rate;

    % 3. Calculate Dynamic Tau
    % tau is inversely proportional to activity.
    % High activity = Small tau (fast decay)
    current_tau = params.K / (state.activity + 1e-5);

    % 4. Apply GLOBAL DECAY to the existing surface
    % This is the mathematically correct step my previous code skipped.
    % We decay the WHOLE surface by the amount of time passed in this batch,
    % using the CURRENT tau.
    decay_factor = exp(-dt_batch / current_tau);
    
    % Apply decay to the background (everything that happened before)
    S = S * decay_factor;

    % 5. Add NEW Events
    % For a standard Time Surface, we set the pixel to 1.0 (Reset).
    % Since we are processing a batch, we need to handle indices carefully.
    
    % Linear indices of new events
    linear_idx = sub2ind([H, W], x_list, y_list);
    
    % Set these pixels to max value (Reset)
    S(linear_idx) = 1.0; 
    
    % 6. Update State
    state.last_update_time = t_batch_end;

    % Normalize the surface
    normalized_output_frame = (S + 1) / 2;
    
end