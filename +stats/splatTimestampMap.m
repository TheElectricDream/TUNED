function [splatted_map] = splatTimestampMap(x, y, t, imgSz, sigma)
% SPLATTIMESTAMPMAP Efficiently spreads event timestamps spatially using 
% Gaussian splatting (Normalized Convolution).
%
% INPUTS:
%   x, y   - Vectors of x and y coordinates (matching image dimensions)
%   t      - Vector of timestamp values
%   imgSz  - [SizeX, SizeY] (e.g., [640, 480])
%   sigma  - Standard deviation of Gaussian kernel (controls "X by X" spread)
%
% OUTPUT:
%   splatted_map - 2D matrix with spatially spread timestamp values.

    % 1. VALIDATION & PRE-ALLOCATION
    % Ensure inputs are double for precision during math, or single for speed
    % Using single saves memory and is faster for GPU/CPU calc.
    if ~isa(t, 'single'), t = single(t); end
    
    % Initialize grids
    raw_map = zeros(imgSz, 'single');
    mask_map = zeros(imgSz, 'single');
    
    % 2. FAST MAPPING (Vectorized)
    % Convert (x,y) to linear indices for fast assignment.
    % Note: MATLAB uses (row, col). If your x is horizontal (1..640) and 
    % y is vertical (1..480), and imgSz is [640, 480], ensure mapping is correct.
    % Based on your main.m: sub2ind(imgSz, x, y) implies 
    % dim 1 is X, dim 2 is Y.
    linear_idx = sub2ind(imgSz, x, y);
    
    % Assign values. 
    % Since 't' is sorted, the later values in the vector will overwrite 
    % earlier ones at the same pixel. This automatically handles "collisions"
    % by keeping the most recent timestamp (which is usually desired).
    raw_map(linear_idx) = t;
    mask_map(linear_idx) = 1.0;
    
    % 3. GAUSSIAN SPLATTING (Normalized Convolution)
    % We use 'Padding', 0 to ensure we don't bleed edge artifacts, 
    % and FilterSize to limit computation to relevant neighborhood.
    
    % Calculate filter size based on sigma to capture ~99% of the curve
    k_size = 2 * ceil(2 * sigma) + 1; 
    
    % Blur the raw values
    blurred_vals = imgaussfilt(raw_map, sigma, ...
        'FilterSize', k_size, 'Padding', 0);
        
    % Blur the mask (weights)
    blurred_weights = imgaussfilt(mask_map, sigma, ...
        'FilterSize', k_size, 'Padding', 0);
    
    % 4. NORMALIZE
    % Divide to recover magnitude. Add eps to avoid divide-by-zero.
    splatted_map = blurred_vals ./ (blurred_weights + eps('single'));
    
    % Optional: Clean up areas that had virtually no events
    % (Floating point division can leave tiny artifacts in empty space)
    splatted_map(blurred_weights < 1e-4) = 0;

end