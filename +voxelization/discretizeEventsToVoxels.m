function [voxelOccupancy, xEdges, yEdges,...
    tEdges] = discretizeEventsToVoxels(x, y, t, opts)
    
    % Define bin edges (voxel boundaries)
    %xEdges = (floor(min(x)/opts.dx)*opts.dx) : opts.dx : (ceil(max(x)/opts.dx)*opts.dx + opts.dx);
    %yEdges = (floor(min(y)/opts.dy)*opts.dy) : opts.dy : (ceil(max(y)/opts.dy)*opts.dy + opts.dy);
    tEdges = (floor(min(t)/opts.dt)*opts.dt) : opts.dt : (ceil(max(t)/opts.dt)*opts.dt + opts.dt);
    
    xEdges = 1:641;
    yEdges = 1:481;


    % Bin each point into voxel indices
    ix = discretize(x, xEdges);
    iy = discretize(y, yEdges);
    it = discretize(t, tEdges);
    
    keep = ~isnan(ix) & ~isnan(iy) & ~isnan(it);
    subs = [ix(keep), iy(keep), it(keep)];
    
    % Build a voxel volume: counts per voxel (dense 3D array)
    sz = [numel(xEdges)-1, numel(yEdges)-1, numel(tEdges)-1];
    Vcount = accumarray(subs, 1, sz, @sum, 0);   % count of points in each voxel
    
    % Use occupancy instead of counts
    voxelOccupancy = Vcount > 0;
    
end