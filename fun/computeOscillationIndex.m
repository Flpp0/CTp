function OI = computeOscillationIndex(r_pad)
% computeOscillationIndex Computes the oscillation index (OI) for a residue function.
%
% Syntax:
%   OI = computeOscillationIndex(r_pad)
%
% Description:
%   This function calculates the oscillation index (OI) as described in Wu et al. (2003). 
%   The oscillation index quantifies the smoothness of the deconvolved residue function, 
%   and is defined as:
%
%       OI = (1/(L * fmax)) * sum(|r(t) - 2*r(t-1) + r(t-2)|) for t = 3,...,L,
%
%   where r_pad is an L-by-1 vector representing the padded residue function and fmax is 
%   the maximum value of r_pad. A higher OI indicates a more oscillatory (less smooth) signal.
%
% Inputs:
%   r_pad - L-by-1 vector of the padded residue function.
%
% Outputs:
%   OI    - Scalar value representing the oscillation index.
%
% Reference:
%   O. Wu et al., “Tracer arrival timing‐insensitive technique for estimating flow in MR 
%   perfusion‐weighted imaging using singular value decomposition with a block‐circulant 
%   deconvolution matrix,” Magnetic Resonance in Med, vol. 50, no. 1, pp. 164–174, Jul. 2003, 
%   doi: 10.1002/mrm.10522.
%
% Example:
%   r_pad = rand(100,1); % Example residue function (padded)
%   OI = computeOscillationIndex(r_pad);
%

    fmax = max(r_pad);
    if fmax <= 0
        OI = 0;
        return;
    end
    L = length(r_pad);
    sumDiff = 0;
    for t = 3:L
        sumDiff = sumDiff + abs(r_pad(t) - 2*r_pad(t-1) + r_pad(t-2));
    end
    OI = sumDiff / (fmax * L);
end
