function [all_pts] = buildSpatioTemporalMap(frame_storage_cells)

    % Build a spatiotemporal point cloud from alts_frame_storage
    all_pts = [];
    for frameIndex = 1:length(frame_storage_cells)
        frame = frame_storage_cells{frameIndex};
        frame(frame==0.5)=nan;
        mask  = single(abs(frame) > 0);
        [rx, ry] = find(mask);
        z = repmat(frameIndex, numel(rx), 1);
        all_pts = [all_pts; rx, ry, z]; %#ok<AGROW>
    end

end