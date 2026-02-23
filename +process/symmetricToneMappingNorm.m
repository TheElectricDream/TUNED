function norm_S = symmetricToneMappingNorm(S, scale)
    % Fixed symmetric tone mapping via hyperbolic tangent.
    % Maps S -> [0, 1] with midpoint at 0.5 (zero surface = gray).
    %
    %   scale controls contrast:
    %     - Larger scale  = softer curve, more headroom for extremes
    %     - Smaller scale = steeper curve, more contrast in quiet regions
    %
    % The tanh function is a standard sigmoidal tone-mapping operator;
    % see Reinhard et al. (2002), "Photographic Tone Reproduction for
    % Digital Images," ACM SIGGRAPH, Eq. 4 and discussion of sigmoid
    % compression for high dynamic range imagery.

    if nargin < 2
        scale = 3.0;
    end

    norm_S = 0.5 + 0.5 * tanh(S / scale);
end
