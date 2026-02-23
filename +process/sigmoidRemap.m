function x_out = sigmoidRemap(x, tau_min, tau_max, k)
%SIGMOID_REMAP Remap values using a sigmoid between specified bounds.
%   X_OUT = SIGMOIDREMAP(X, TAU_MIN, TAU_MAX) remaps the input array X to
%   the range [TAU_MIN, TAU_MAX] using a sigmoid-shaped mapping. The input
%   values are first normalized to [-k, k] (default k=6) before applying
%   the sigmoid.
%
%   X_OUT = SIGMOID_REMAP(X, TAU_MIN, TAU_MAX, K) specifies the steepness
%   parameter K used for normalization. Larger K produces a steeper
%   transition in the sigmoid.
%
%   Inputs:
%     X        - Numeric vector or array of values to remap.
%     TAU_MIN  - Scalar specifying the lower bound of the output range.
%     TAU_MAX  - Scalar specifying the upper bound of the output range.
%     K        - (Optional) Scalar positive parameter controlling the
%                sigmoid steepness. Default is 6.
%
%   Output:
%     X_OUT    - Numeric array of same size as X with values mapped to the
%                interval [TAU_MIN, TAU_MAX].
%
%   Example:
%     y = sigmoid_remap(randn(100,1), 0, 1, 8);
%
%   Notes:
%     - If X is constant, the function returns a constant array equal to
%       the midpoint (TAU_MIN + TAU_MAX)/2 to avoid division by zero.
%
%   See also EXP.

    if nargin < 4
        k = 6;
    end
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
    x_norm = (x - xmin) / (xmax - xmin) * 2 * k - k;
    s = 1 ./ (1 + exp(-x_norm));
    x_out = tau_min + s .* (tau_max - tau_min);
end