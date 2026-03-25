function [S, state, normalized_output_frame] = ...
    adaptiveGlobalDecay(S, x_list, y_list, t_list, state, params)
% ADAPTIVEGLOBALDECAY  Adaptive Global Decay (AGD) time surface.
%
%   [S, STATE, NORMALIZED_OUTPUT_FRAME] = ADAPTIVEGLOBALDECAY(S,
%   X_LIST, Y_LIST, T_LIST, STATE, PARAMS) updates the surface using
%   a globally adaptive decay rate derived from the instantaneous
%   event rate. High event activity produces fast decay (small tau);
%   low activity produces slow decay (long memory).
%
%   Reference:
%       Nunes et al., "Adaptive Global Decay Process for Event
%       Cameras," Proc. IEEE CVPR, pp. 9771-9780, 2023.
%       DOI: 10.1109/CVPR52729.2023.00942
%
%   Inputs:
%     S      - [imgSz] Current surface state.
%     x_list - [N x 1] Row coordinates of events.
%     y_list - [N x 1] Column coordinates of events.
%     t_list - [N x 1] Timestamps [s], sorted ascending.
%     state  - Struct with persistent fields:
%                .activity         - Smoothed event rate (scalar)
%                .last_update_time - Timestamp of previous batch
%     params - Struct with tuning parameters:
%                .alpha - EMA smoothing factor for activity
%                .K     - Scaling constant: tau = K / activity
%
%   Outputs:
%     S                       - [imgSz] Updated surface.
%     state                   - Updated state struct.
%     normalized_output_frame - [imgSz] Display surface in [0, 1].
%
%   Algorithm:
%     1. Compute global event rate for this batch.
%     2. Smooth the rate with a leaky integrator (EMA).
%     3. Derive dynamic tau = K / activity.
%     4. Decay the entire surface by exp(-dt / tau).
%     5. Reset pixels at event locations to 1.0.
%
%   Notes:
%     - Reset-based: event pixels are set to 1.0 (overwrite).
%     - tau is global — all pixels share the same decay rate.
%     - Unsigned output in [0, 1]. Direct passthrough.
%     - Coordinates: x = row, y = col, sub2ind(imgSz, x, y).
%
%   See also: accumulator.timeSurface,
%             accumulator.localAdaptiveTimeSurface

    % ----------------------------------------------------------------
    % 0. Early exit
    % ----------------------------------------------------------------
    if isempty(t_list)
        normalized_output_frame = S;
        return;
    end

    % ----------------------------------------------------------------
    % 1. Compute batch time delta
    % ----------------------------------------------------------------
    t_batch_end = t_list(end);
    dt_batch = t_batch_end - state.last_update_time;

    % Guard against zero or negative dt (duplicate/out-of-order)
    if dt_batch < 1e-6
        fprintf('\ndt_batch is zero, setting to 1e-6\n');
    end
    dt_batch = max(dt_batch, 1e-6);

    % ----------------------------------------------------------------
    % 2. Estimate global event activity (leaky integrator)
    % ----------------------------------------------------------------
    current_rate = numel(t_list) / dt_batch;
    state.activity = (1 - params.alpha) * state.activity ...
                   + params.alpha * current_rate;

    % ----------------------------------------------------------------
    % 3. Compute dynamic tau (inversely proportional to activity)
    % ----------------------------------------------------------------
    current_tau = params.K / (state.activity + 1e-5);

    % ----------------------------------------------------------------
    % 4. Apply global decay to the entire surface
    % ----------------------------------------------------------------
    decay_factor = exp(-dt_batch / current_tau);
    S = S * decay_factor;

    % ----------------------------------------------------------------
    % 5. Reset event pixels to maximum value
    % ----------------------------------------------------------------
    [H, W] = size(S);
    linear_idx = sub2ind([H, W], x_list, y_list);
    S(linear_idx) = 1.0;

    % ----------------------------------------------------------------
    % 6. Update state and normalize for display
    % ----------------------------------------------------------------
    state.last_update_time = t_batch_end;

    % Surface is unsigned [0, 1] — output directly.
    normalized_output_frame = S;

end