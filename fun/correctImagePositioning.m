function registeredHUImages = correctImagePositioning(smoothedHUImages, registeredHUImages, uniqueTimePointsSeconds, upsample_factor, highMovementIndices)
% correctImagePositioning Adjusts the positioning of registered images to match original slice order.
%
% Syntax:
%   registeredHUImages = correctImagePositioning(smoothedHUImages, registeredHUImages, uniqueTimePointsSeconds, upsample_factor)
%   registeredHUImages = correctImagePositioning(smoothedHUImages, registeredHUImages, uniqueTimePointsSeconds, upsample_factor, highMovementIndices)
%
% Description:
%   This function reassigns slices from the registered (upsampled) volumes so that their order corresponds
%   to the original (smoothed) images. A set of valid volumes is computed by excluding time points
%   with high movement (if provided).
%
% Inputs:
%   smoothedHUImages      - Cell array of original images (nSlices x nTimePoints).
%   registeredHUImages    - Cell array of registered (upsampled) images (nUpsampledSlices x nTimePoints).
%   uniqueTimePointsSeconds - Numeric array of time points (in seconds) corresponding to the original data.
%   upsample_factor       - Upsampling factor (1 means no upsampling).
%   highMovementIndices   - (Optional) Indices to exclude (e.g., high movement); default is [].
%
% Output:
%   registeredHUImages    - Corrected cell array of registered images with slice ordering matching the original.
%
% Example:
%   regImages = correctImagePositioning(smoothedHUImages, registeredHUImages, uniqueTimePointsSeconds, 1);
%   regImages = correctImagePositioning(smoothedHUImages, registeredHUImages, uniqueTimePointsSeconds, 1, highMovementIndices);


    if nargin < 5
        highMovementIndices = [];
    end

    numOriginalSlices = size(smoothedHUImages, 1);
    numUpsampledSlices = size(registeredHUImages, 1);
    registeredHUImagesSorted = cell(size(smoothedHUImages));

    numTimePointsOriginal = size(smoothedHUImages, 2);
    validVolumes = setdiff(1:numTimePointsOriginal, highMovementIndices);
    
    for sIdx = 1:numOriginalSlices
        slice_idx_in_registeredHUImages = round((sIdx - 1) * upsample_factor + 1);
        if slice_idx_in_registeredHUImages > numUpsampledSlices
            warning('Slice index %d exceeds number of upsampled slices (%d).', slice_idx_in_registeredHUImages, numUpsampledSlices);
            continue;
        end
        
        currentSIdxTimes = find(~cellfun('isempty', smoothedHUImages(sIdx, :)));
        validIdx = validVolumes(validVolumes <= length(currentSIdxTimes));
        if isempty(validIdx)
            continue;
        end
        
        registeredHUImagesSorted(sIdx, currentSIdxTimes(validIdx)) = registeredHUImages(slice_idx_in_registeredHUImages, validIdx);
    end
    
    registeredHUImages = registeredHUImagesSorted;
end
