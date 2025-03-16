function [R_all] = computeResidueFunctions_cSVD(AIF, TCCs, lambda, m, dt)
% computeResidueFunctions_cSVD Computes flow-scaled residue functions using block-circulant cSVD.
%
% Syntax:
%   R_all = computeResidueFunctions_cSVD(AIF, TCCs, lambda, m, dt)
%
% Description:
%   This function performs deconvolution for multiple voxels using a block-circulant 
%   matrix representation of the convolution process. After applying Simpson-like filtering 
%   to the AIF, it constructs a circulant matrix and computes its SVD to obtain a pseudo-inverse.
%   The Tissue Attenuation Curves (TACs) are zero-padded, deconvolved, scaled by 1/dt, and then
%   cropped back to the original time length.
%
%   This method is based on the technique described in:
%
%       O. Wu et al., “Tracer arrival timing‐insensitive technique for estimating flow in MR 
%       perfusion‐weighted imaging using singular value decomposition with a block‐circulant 
%       deconvolution matrix,” Magnetic Resonance in Med, vol. 50, no. 1, pp. 164–174, Jul. 2003, 
%       doi: 10.1002/mrm.10522.
%
% Inputs:
%   AIF    - 1 x T vector representing the arterial input function.
%   TCCs   - V x T matrix of tissue concentration curves (V = number of voxels).
%   lambda - Scalar; fraction of the largest singular value for thresholding.
%   m      - Extension factor for the block-circulant matrix.
%   dt     - Time sampling interval (in seconds).
%
% Output:
%   R_all  - V x T matrix of residue functions scaled by cerebral blood flow.
%
% See also: gallery, svd


    % Step 1: Simpson-like filtering on the AIF
    T = size(AIF,2);
    AIF_filt = AIF;
    if T > 2
        for k = 2:(T-1)
            AIF_filt(k) = (AIF(k-1) + 4*AIF(k) + AIF(k+1)) / 6;
        end
    end
    
    % Step 2: Zero-pad the AIF and build the block-circulant matrix
    L = m * T;
    Ca_pad = zeros(1, L);
    Ca_pad(1:T) = AIF_filt;
    D = gallery('circul', Ca_pad)';  % L x L
    
    % Step 3: Perform SVD on the circulant matrix and truncate
    [U, S, Vmat] = svd(D);
    sigmaMax = max(diag(S));
    keepMask = diag(S) >= (lambda * sigmaMax);
    invS = diag(1 ./ diag(S));
    invS(~keepMask, ~keepMask) = 0;
    
    % Step 4: Zero-pad TCCs
    Vvox = size(TCCs,1);
    C_pad = zeros(L, Vvox);
    C_pad(1:T, :) = TCCs';  % shape: [L x Vvox]
    
    % Step 5: Deconvolve in one shot, then multiply by (1/dt)
    R_all_pad = (1/dt) * (Vmat * invS * U' * C_pad);
    
    % Step 6: Crop back to original length
    R_cropped = R_all_pad(1:T, :);  % [T x Vvox]
    
    % Step 7: Return as V x T
    R_all = R_cropped';
end
