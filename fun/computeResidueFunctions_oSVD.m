function [R_all, OI_all] = computeResidueFunctions_oSVD(AIF, TCCs, dt, OI_limit, m, alphaCandidates)
% computeResidueFunctions_oSVD Computes residue functions using oscillation-index SVD (oSVD).
%
% Syntax:
%   [R_all, OI_all] = computeResidueFunctions_oSVD(AIF, TCCs, dt, OI_limit, m, alphaCandidates)
%
% Description:
%   This function implements an oscillation-index (OI) based deconvolution using a block-circulant
%   matrix representation of the convolution operator. The AIF is first filtered using Simpson-like filtering,
%   and a circulant matrix is constructed. For each voxel, a range of alpha thresholds is explored 
%   (given by alphaCandidates) to compute candidate residue functions. The first candidate for which the
%   oscillation index is below OI_limit is locked in; if none are below the threshold, the candidate with the
%   lowest OI is selected. Finally, the solution is cropped back to the original time length.

%   This method is based on the technique described in:
%
%       O. Wu et al., “Tracer arrival timing‐insensitive technique for estimating flow in MR 
%       perfusion‐weighted imaging using singular value decomposition with a block‐circulant 
%       deconvolution matrix,” Magnetic Resonance in Med, vol. 50, no. 1, pp. 164–174, Jul. 2003, 
%       doi: 10.1002/mrm.10522.
%
% Inputs:
%   AIF             - 1 x T vector representing the arterial input function.
%   TCCs            - V x T matrix of tissue concentration curves.
%   dt              - Time sampling interval (seconds).
%   OI_limit        - Maximum acceptable oscillation index.
%   m               - Extension factor for the block-circulant matrix.
%   alphaCandidates - Array of candidate alpha thresholds (e.g., 0.01:0.01:0.3).
%
% Outputs:
%   R_all  - V x T matrix of flow-scaled residue functions computed via oSVD.
%   OI_all - V x 1 vector of final oscillation index values for each voxel.
%
% See also: gallery, svd, computeOscillationIndex

    % Step 1: Simpson-like filtering
    T = length(AIF);
    AIF_filt = AIF;
    if T > 2
        for k = 2:(T-1)
            AIF_filt(k) = (AIF(k-1) + 4*AIF(k) + AIF(k+1)) / 6;
        end
    end
    
    % Step 2: Build block-circulant matrix with zero-padding
    L = m * T;
    Ca_pad = zeros(1, L);
    Ca_pad(1:T) = AIF_filt;
    D = gallery('circul', Ca_pad)';
    
    % Step 3: SVD of the circulant matrix
    [U,S,V] = svd(D);
    sigmaVec = diag(S);
    sigmaMax = max(sigmaVec);
    invS_full = diag(1 ./ sigmaVec);
    
    % Step 4: Zero-pad TCCs
    Vvox = size(TCCs,1);
    C_pad = zeros(L, Vvox);
    C_pad(1:T,:) = TCCs.';  % [L x Vvox]
    
    % Step 5: Determine optimal alpha for each voxel via oscillation index
    nAlpha = numel(alphaCandidates);
    locked   = false(Vvox,1);
    bestOI   = inf(Vvox,1);
    bestRpad = zeros(L, Vvox);
    
    for aIdx = 1:nAlpha
        alpha = alphaCandidates(aIdx);
        keepMask = (sigmaVec >= alpha * sigmaMax);
        W = invS_full;
        W(~keepMask, ~keepMask) = 0;
        invD_local = V * W * U';
        R_candidate = invD_local * C_pad;
        R_candidate = (1/dt) * R_candidate;
        
        for v = 1:Vvox
            if ~locked(v)
                r_voxel = R_candidate(:,v);
                OI_val  = computeOscillationIndex(r_voxel);
                if OI_val < bestOI(v)
                    bestOI(v) = OI_val;
                    bestRpad(:,v) = r_voxel;
                end
                if OI_val < OI_limit
                    locked(v) = true;
                end
            end
        end
        
        if all(locked)
            break;
        end
    end
    
    R_all_pad = bestRpad;
    OI_all = bestOI;
    R_cropped = R_all_pad(1:T, :);
    R_all = R_cropped.';
end
