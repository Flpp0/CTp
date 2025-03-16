function [CBF_map, CBV_map, MTT_map, TMAX_map] = computePerfusionMaps(residueAll, tissueAll, time, mask, AIF)
% computePerfusionMaps Computes voxel-wise perfusion parameters (CBF, CBV, MTT, TMAX) from deconvolved tissue data.
%
% Syntax:
%   [CBF_map, CBV_map, MTT_map, TMAX_map] = computePerfusionMaps(residueAll, tissueAll, time, mask, AIF)
%
% Description:
%   This function calculates perfusion maps from the deconvolved residue functions and the raw
%   tissue attenuation curves. Given:
%     - residueAll: a [V x T] matrix containing the deconvolved residue function r(t) for each voxel,
%     - tissueAll:  a [V x T] matrix containing the raw tissue attenuation curves,
%     - time:       a 1 x T vector of time points (in seconds),
%     - mask:       a logical array indicating brain voxels, and
%     - AIF:        a 1 x T vector representing the arterial input function,
%
%   the function computes:
%
%       CBV = ( ∫ tissue(t) dt / ∫ AIF(t) dt ) * Hct * (1/rho) * 100       [ml/100g]
%       CBF = max(r(t)) * Hct * (1/rho) * 60 * 100                         [ml/100g/min]
%       MTT = (CBV / CBF) * 60                                             [s]
%       TMAX = argmax(r(t))                                                [s]
%
%   where Hct (hematocrit) is approximately 0.73, rho (tissue density) is 1.04 g/mL, and time is in seconds.
%
% Inputs:
%   residueAll - [V x T] matrix of deconvolved residue functions.
%   tissueAll  - [V x T] matrix of raw tissue attenuation curves.
%   time       - 1 x T vector of time points (seconds).
%   mask       - Logical array indicating brain voxels.
%   AIF        - 1 x T arterial input function.
%
% Outputs:
%   CBF_map   - Perfusion map for Cerebral Blood Flow (ml/100g/min), same dimensions as mask.
%   CBV_map   - Perfusion map for Cerebral Blood Volume (ml/100g), same dimensions as mask.
%   MTT_map   - Perfusion map for Mean Transit Time (s), same dimensions as mask.
%   TMAX_map  - Map indicating the time (in seconds) at which r(t) is maximal, same dimensions as mask.
%
% Example:
%   [CBF_map, CBV_map, MTT_map, TMAX_map] = computePerfusionMaps(R_all, TCC_all, time, mask, AIF);
%
% Reference:
%   A. Fieselmann, M. Kowarschik, A. Ganguly, J. Hornegger, and R. Fahrig, “Deconvolution-Based
%   CT and MR Brain Perfusion Measurement: Theoretical Model Revisited and Practical Implementation
%   Details,” International Journal of Biomedical Imaging, vol. 2011, pp. 1–20, 2011, 
%   doi: 10.1155/2011/467563.


    if ~islogical(mask)
        error('mask must be a logical array.');
    end

    [numVox_r,  numTime_r ] = size(residueAll);
    [numVox_tc, numTime_tc] = size(tissueAll);

    if (numVox_r~=numVox_tc) || (numTime_r~=numTime_tc)
        error('Dimension mismatch: residueAll vs. tissueAll');
    end
     nVox = numVox_r;

    if sum(mask(:)) ~= nVox
        error('mask dimension mismatch with residueAll/tissueAll');
    end
    if length(time) ~= numTime_r
        error('Time vector does not match # of columns in residueAll/tissueAll');
    end
    if length(AIF) ~= numTime_r
        error('AIF length does not match # of time points in residueAll/tissueAll');
    end

    % Constants
    Hct = 0.73;   
    rho = 1.04;   % g/mL
    intAIF = trapz(time, AIF);

    % Preallocate results
    CBF_map   = zeros(size(mask));
    CBV_map   = zeros(size(mask));
    MTT_map   = zeros(size(mask));
    TMAX_map  = zeros(size(mask));

    brainVox = find(mask);  % linear indices of mask

    for i = 1:length(brainVox)
        voxIdx  = brainVox(i);
        r_vox   = residueAll(i, :);   % 1D
        c_vox   = tissueAll(i, :);    % 1D

        % 1) Find r_max and tMax
        [r_max, idxMax] = max(r_vox);
        tMax_sec = time(idxMax);   % TMAX in seconds

        % 2) Convert r_max => CBF => ml/100g/min
        %    Multiply by: (1/rho)*60*100
        CBF = r_max * Hct * (1/rho) * 60 * 100;

        % 3) Tissue integral for CBV
        intTissue = trapz(time, c_vox);
        if intAIF > 0
            cbv_ratio = intTissue / intAIF;
        else
            cbv_ratio = 0;
        end

        % 4) CBV => ml/100g
        %    multiply by Hct*(1/rho)*100
        CBV = cbv_ratio * Hct * (1/rho) * 100;

        % 5) MTT => seconds
        %    MTT = (CBV / CBF)*60 => if CBF is in ml/100g/min
        if CBF>0
            MTT_sec = (CBV / CBF) * 60; 
        else
            MTT_sec = 0;
        end

        % Store in maps
        CBF_map(voxIdx)   = CBF;
        CBV_map(voxIdx)   = CBV;
        MTT_map(voxIdx)   = MTT_sec;
        TMAX_map(voxIdx)  = tMax_sec;
    end
end
