function [R_all] = computeResidueFunctions_sSVD(AIF, TCCs, lambdaGlobal, dt)
% computeResidueFunctions_sSVD Computes flow-scaled residue functions using simple truncated SVD.
%
% Syntax:
%   R_all = computeResidueFunctions_sSVD(AIF, TCCs, lambdaGlobal, dt)
%
% Description:
%   This function performs deconvolution on Tissue Attenuation Curves (TACs) using a standard
%   linear convolution matrix constructed from a filtered arterial input function (AIF). A singular
%   value decomposition (SVD) is applied to compute the pseudo-inverse, after which a truncation is
%   performed based on a fraction of the largest singular value. The resulting residue functions are
%   then scaled by 1/dt.
%
%   This approach is primarily based on:
%
%       A. Fieselmann, M. Kowarschik, A. Ganguly, J. Hornegger, and R. Fahrig, “Deconvolution-Based
%       CT and MR Brain Perfusion Measurement: Theoretical Model Revisited and Practical Implementation
%       Details,” International Journal of Biomedical Imaging, vol. 2011, pp. 1–20, 2011, 
%       doi: 10.1155/2011/467563.
%
%       O. Wu et al., “Tracer arrival timing‐insensitive technique for estimating flow in MR 
%       perfusion‐weighted imaging using singular value decomposition with a block‐circulant 
%       deconvolution matrix,” Magnetic Resonance in Med, vol. 50, no. 1, pp. 164–174, Jul. 2003, 
%       doi: 10.1002/mrm.10522.
%
% Inputs:
%   AIF          - 1 x T vector representing the arterial input function.
%   TCCs         - V x T matrix of tissue concentration curves.
%   lambdaGlobal - Scalar; fraction of the largest singular value for thresholding.
%   dt           - Time sampling interval (in seconds).
%
% Output:
%   R_all        - V x T matrix of flow-scaled residue functions.
%
% See also: svd


    [~, T] = size(AIF);
    
    % Step 1: Simpson-like filtering
    AIF_filt = AIF;
    if T > 2
        for k = 2:(T-1)
            AIF_filt(k) = (AIF_filt(k-1) + 4*AIF_filt(k) + AIF_filt(k+1)) / 6;
        end
    end
    
    % Step 2: Build standard linear convolution matrix A (T x T).
    A = zeros(T, T);
    for i = 1:T
        for j = 1:i
            A(i,j) = AIF_filt(i-j+1);
        end
    end
    
    % Step 3: SVD of A
    [U, S, Vmat] = svd(A);
    sigmaVec = diag(S);
    sigmaMax = max(sigmaVec);
    
    % Step 4: Truncate
    keepMask = sigmaVec >= (lambdaGlobal * sigmaMax);
    invS = diag(1 ./ sigmaVec);
    invS(~keepMask, ~keepMask) = 0;
    
    % Step 5: Build pseudo-inverse
    invA = Vmat * invS * U';
    
    % Step 6: Prepare TCCs for multiplication
    C_mat = TCCs';  % shape (T x V)
    
    % R_candidate = invA * C_mat => (T x V)
    R_candidate = invA * C_mat;
    % Then scale by (1/dt) if you want residue in 1/s
    R_candidate = (1/dt) * R_candidate;
    
    % Step 7: Return as (V x T)
    R_all = R_candidate';
end
