function x_out = linearRemap(x, tau_min, tau_max)
%LINEARREMAP Remap values linearly between specified bounds.
%   X_OUT = LINEARREMAP(X, TAU_MIN, TAU_MAX) remaps the input array X to
%   the range [TAU_MIN, TAU_MAX] using a linear mapping.
%
%   Inputs:
%     X        - Numeric vector or array of values to remap.
%     TAU_MIN  - Scalar specifying the lower bound of the output range.
%     TAU_MAX  - Scalar specifying the upper bound of the output range.
%
%   Output:
%     X_OUT    - Numeric array of same size as X with values mapped to the
%                interval [TAU_MIN, TAU_MAX].
%
%   Example:
%     y = linearRemap(randn(100,1), 0, 1);
%
%   Notes:
%     - If X is constant, the function returns a constant array equal to
%       the midpoint (TAU_MIN + TAU_MAX)/2 to avoid division by zero.
%
%   See also SIGMOIDREMAP.

    x = double(x);

    if isempty(x)
        x_out = zeros(size(x));
        return;
    end

    xmin = min(x(:));
    xmax = max(x(:));

    if xmax == xmin
        x_out = repmat((tau_min + tau_max) / 2, size(x));
        return;
    end

    x_norm = (x - xmin) / (xmax - xmin);
    x_out = tau_min + x_norm .* (tau_max - tau_min);
end