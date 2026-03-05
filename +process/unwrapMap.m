function [filtered_events] = unwrapMap(...
    map, sorted_x, sorted_y, sorted_t, t_offset) %#ok<*DEFNU>
% UNWRAPMAP  Extract events passing a coherence threshold from a 2D map.
%
%   FILTERED_EVENTS = UNWRAPMAP(MAP, SORTED_X, SORTED_Y, SORTED_T,
%   T_OFFSET) looks up the coherence score at each event's pixel
%   location and returns only those events where the score is
%   nonzero (i.e., passing the upstream threshold).
%
%   Inputs:
%     map      - [imgSz] 2D coherence score map (pre-thresholded).
%     sorted_x - [N x 1] Row coordinates of events in current frame.
%     sorted_y - [N x 1] Column coordinates of events.
%     sorted_t - [N x 1] Timestamps [s] relative to frame start.
%     t_offset - Scalar time offset to convert to absolute time.
%
%   Outputs:
%     filtered_events - [M x 4] Array of [x, y, t_abs, score] for
%                       events that pass the coherence threshold.
%
%   Notes:
%     - Pre-allocates output for performance (avoids growing array).
%     - Coordinates: x = row, y = col. The map is indexed as
%       map(x, y) following the repository convention.
%
%   See also: coherence.computeCoherenceMask

    % ----------------------------------------------------------------
    % 0. Pre-allocate output (worst case: all events pass)
    % ----------------------------------------------------------------
    n_events = length(sorted_x);
    filtered_events = zeros(n_events, 4);
    count = 0;

    % ----------------------------------------------------------------
    % 1. Look up coherence score for each event
    % ----------------------------------------------------------------
    for i = 1:n_events
        x = sorted_x(i);
        y = sorted_y(i);
        t = sorted_t(i);

        coherence_score = map(x, y);

        if coherence_score > 0
            count = count + 1;
            filtered_events(count, :) = ...
                [x, y, t + t_offset, coherence_score];
        end
    end

    % Trim unused rows
    filtered_events = filtered_events(1:count, :);

end