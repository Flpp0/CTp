function [uniqueTimePoints, uniqueTimePointsSeconds, uniqueSliceLocations, sortedHUImages, sortedPixelSpacing] = sortImages(vpctImages, timePoints, sliceLocations, vpctPixelSpacing)
% sortImages Sorts VPCT images into a cell array organized by slice and time.
%
% Syntax:
%   [uniqueTimePoints, uniqueTimePointsSeconds, uniqueSliceLocations, sortedHUImages, sortedPixelSpacing] = sortImages(vpctImages, timePoints, sliceLocations, vpctPixelSpacing)
%
% Description:
%   This function identifies unique time points (converting them to seconds) and slice locations.
%   It then arranges the VPCT images into a 2D cell array where rows correspond to slices and columns
%   correspond to unique time points.
%
% Inputs:
%   vpctImages     - Cell array of VPCT images.
%   timePoints     - Cell array of acquisition times.
%   sliceLocations - Numeric array of slice locations.
%   vpctPixelSpacing- Cell array of PixelSpacing values.
%
% Outputs:
%   uniqueTimePoints         - Cell array of unique acquisition times.
%   uniqueTimePointsSeconds  - Numeric array of unique times in seconds.
%   uniqueSliceLocations     - Numeric array of unique slice locations.
%   sortedHUImages           - 2D cell array (nSlices x nTimePoints) of images.
%   sortedPixelSpacing       - Cell array of PixelSpacing values per slice.
%
% Example:
%   [uniqueTimePoints, uniqueTimePointsSeconds, uniqueSliceLocations, sortedHUImages, sortedPixelSpacing] = ...
%         sortImages(vpctImages, timePoints, sliceLocations, vpctPixelSpacing);
%
% See also: timeStringToSeconds

    uniqueTimePoints = unique(timePoints);
    uniqueTimePointsSeconds = zeros(size(uniqueTimePoints));
    for time_idx = 1:length(uniqueTimePointsSeconds)
        uniqueTimePointsSeconds(time_idx) = timeStringToSeconds(uniqueTimePoints{time_idx});
    end
    uniqueSliceLocations = unique(sliceLocations);
    numSlices = length(uniqueSliceLocations);
    sortedHUImages = cell(numSlices, length(uniqueTimePoints));
    sortedPixelSpacing = cell(numSlices, 1);
    
    for i = 1:length(vpctImages)
        tIdx = find(strcmp(uniqueTimePoints, timePoints{i}));
        sIdx = find(uniqueSliceLocations == sliceLocations(i));
        sortedHUImages{sIdx, tIdx} = vpctImages{i};
        if isempty(sortedPixelSpacing{sIdx})
            sortedPixelSpacing{sIdx} = vpctPixelSpacing{i};
        end
    end
end
