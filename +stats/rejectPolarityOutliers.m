function [map, map_raw] = rejectPolarityOutliers(map, map_raw, n_sigma)
% REJECTPOLARITYOUTLIERS  Per-polarity sigma-clipping to median.
%
%   [MAP, MAP_RAW] = REJECTPOLARITYOUTLIERS(MAP, MAP_RAW, N_SIGMA)
%   clamps values beyond N_SIGMA standard deviations from the
%   polarity-conditional mean, replacing them with the polarity
%   median. Both MAP and MAP_RAW are clipped using the same mask
%   (derived from MAP) so that display and feedback paths stay
%   consistent.

    for pol = [1, -1]
        if pol == 1
            mask = map > 0;
        else
            mask = map < 0;
        end

        vals = map(mask);
        if isempty(vals), continue; end

        mu  = mean(vals);
        sig = std(vals);
        med = median(vals);

        if pol == 1
            outliers = map > (mu + n_sigma * sig);
        else
            outliers = map < (mu - n_sigma * sig);
        end

        map(outliers)     = med;
        map_raw(outliers) = median(map_raw(mask));
    end
end
