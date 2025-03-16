function [meanSlices, brainMasks] = performSkullStripping(registeredHUImages, baseDir, patientCode, upsample_factor, nRows, nCols, skullThreshold)
% performSkullStripping Performs skull stripping and computes mean slices.
%
% Syntax:
%   [meanSlices, brainMasks] = performSkullStripping(registeredHUImages, baseDir, patientCode, upsample_factor, nRows, nCols, skullThreshold)
%
% Description:
%   This function reconstructs volumes from registered images, computes mean slices 
%   (pre-contrast images), and performs skull stripping to generate binary brain masks.
%   The skull stripping function is based on the method described in:
%
%       Najm, M., et al. (2019). Automated brain extraction from head CT and CTA images 
%       using convex optimization with shape propagation. Computer Methods and Programs 
%       in Biomedicine, 176, 1–8. doi:10.1016/j.cmpb.2019.04.030.
%
%   At the end of the processing, a montage is created displaying all slices with the brain 
%   mask overlaid in red, and the montage is saved.
%
% Inputs:
%   registeredHUImages - Cell array of registered images (nSlices x nTimePoints).
%   baseDir            - Base directory for patient data.
%   patientCode        - Patient code used in folder names.
%   upsample_factor    - Upsampling factor.
%   nRows, nCols       - Dimensions of the images (e.g., 512, 512).
%   skullThreshold     - Threshold value for skull stripping.
%
% Outputs:
%   meanSlices         - 3D array of mean slices (nRows x nCols x nSlices).
%   brainMasks         - 3D binary mask (nRows x nCols x nSlices) representing the brain region.
%
% Example:
%   [meanSlices, brainMasks] = performSkullStripping(registeredHUImages, baseDir, patientCode, upsample_factor, 512, 512, 100);
%
% See also: SkullStripping (Najm et al., 2019)

    saveDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)]);
    meanSlicesDir = fullfile(saveDir, 'MeanSlices');
    brainMasksDir = fullfile(saveDir, 'BrainMasks');
    meanSlicesFile = fullfile(meanSlicesDir, 'MeanSlices.mat');
    brainMasksFile = fullfile(brainMasksDir, 'BrainMasks.mat');
    
    if ~exist(meanSlicesDir, 'dir'), mkdir(meanSlicesDir); end
    if ~exist(brainMasksDir, 'dir'), mkdir(brainMasksDir); end
    
    numSlices = size(registeredHUImages,1);
    
    if exist(meanSlicesFile, 'file')
        load(meanSlicesFile, 'meanSlices');
        disp("Loading existing meanSlices...");
    else
        disp("Calculating meanSlices for each slice...");
        meanSlices = zeros(nRows, nCols, numSlices);
        numVolumesPreContrast = 4;
        for sIdx = 1:numSlices
            preContrastVolumes = zeros(nRows, nCols, numVolumesPreContrast);
            validTimes = find(~cellfun('isempty', registeredHUImages(sIdx, :)));
            if length(validTimes) < numVolumesPreContrast
                error('Not enough valid pre-contrast volumes for slice %d.', sIdx);
            end
            for tIdx = 1:numVolumesPreContrast
                currentValidTime = validTimes(tIdx);
                preContrastVolumes(:,:,tIdx) = registeredHUImages{sIdx, currentValidTime};
            end
            meanSlice = mean(preContrastVolumes, 3);
            meanSlices(:,:,sIdx) = meanSlice;
        end
        save(meanSlicesFile, 'meanSlices');
        disp("MeanSlices calculated and saved.");
    end
    
    if exist(brainMasksFile, 'file')
        load(brainMasksFile, 'brainMasks');
        disp("Loading existing brain masks...");
    else
        disp("Calculating brain masks using SkullStripping...");
        brainMaskVolume = SkullStripping(double(meanSlices), skullThreshold);
        brainMasks = brainMaskVolume > 0;
        save(brainMasksFile, 'brainMasks');
        disp("Brain masks calculated and saved.");
    end
    
    % Create Montage with All Slices and Brain Mask Overlay
    
    montageDir = fullfile(saveDir, 'SkullStrippingMontage');
    if ~exist(montageDir, 'dir'), mkdir(montageDir); end
    
    overlayImages = zeros(nRows, nCols, 3, numSlices);
    alpha = 0.4; % Transparency factor for overlay
    
    for sIdx = 1:numSlices
        
        sliceImage = mat2gray(meanSlices(:,:,sIdx));
        % Create an RGB image 
        rgbSlice = repmat(sliceImage, [1 1 3]);
        
        currentMask = brainMasks(:,:,sIdx);
        % Using the red channel for the overlay
        overlayImage = rgbSlice; 
        overlayImage(:,:,1) = (1 - alpha) * rgbSlice(:,:,1) + alpha * double(currentMask);
        overlayImage(:,:,2) = (1 - alpha) * rgbSlice(:,:,2);
        overlayImage(:,:,3) = (1 - alpha) * rgbSlice(:,:,3);
        
        overlayImages(:,:,:,sIdx) = overlayImage;
    end
    
    % Figure
    figMontage = figure('Visible','off');
    montage(overlayImages);
    title('Montage of All Slices with Brain Mask Overlay', 'FontSize', 14);
    % Save
    montageFileName = fullfile(montageDir, 'AllSlicesMontage.png');
    exportgraphics(gca, montageFileName, 'Resolution', 300);
    close(figMontage);
end
