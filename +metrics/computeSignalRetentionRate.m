function SRR = computeSignalRetentionRate(n_passed, n_total)
% COMPUTESIGNALRETENTIONRATE  Fraction of events passing coherence filter.
%
%   SRR = COMPUTESIGNALRETENTIONRATE(N_PASSED, N_TOTAL) computes the
%   Signal Retention Rate as the ratio of coherent events to total input
%   events within a processing frame. A guard against division by zero is
%   applied when no events are present.
%
%   Inputs:
%     n_passed - Scalar. Number of events passing the coherence filter.
%     n_total  - Scalar. Total number of input events in the frame.
%
%   Outputs:
%     SRR - Scalar in [0, 1]. Signal retention rate. Returns 0 when
%           n_total is zero.
%
%   Notes:
%     - When n_total == 0, SRR is defined as 0 (no signal to retain).
%     - Both inputs are expected to be non-negative integers, though no
%       type enforcement is applied.
%
%   See also: coherence.computeCoherenceMask

    % ================================================================
    % 0. Handle edge cases
    % ================================================================
    if n_total == 0
        SRR = 999;
        return;
    end

    % ----------------------------------------------------------------
    % 1. Compute retention rate (guarded against empty frames)
    % ----------------------------------------------------------------
    SRR = n_passed / max(n_total, 1);

end