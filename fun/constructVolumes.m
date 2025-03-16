function [volumes, nonEmptyTimeIndices, minTimePoints] = constructVolumes(smoothedHUImages)
% constructVolumes Constructs 3D volumes by aligning non-empty images across slices.
%
% Syntax:
%   [volumes, nonEmptyTimeIndices, minTimePoints] = constructVolumes(smoothedHUImages)
%
% Description:
%   This function examines the cell array of smoothed images (organized by slice and time)
%   and builds a cell array of 3D volumes. Each volume is constructed by taking the nth non-empty
%   image from each slice. It also returns the indices of non-empty time points per slice and the
%   minimum number of non-empty time points.
%
% Inputs:
%   smoothedHUImages - Cell array (nSlices x nTimePoints) of smoothed images.
%
% Outputs:
%   volumes            - Cell array of 3D volumes (each volume is a 3D array).
%   nonEmptyTimeIndices- Cell array containing indices of non-empty images for each slice.
%   minTimePoints      - Minimum number of non-empty images across slices.
%
% Example:
%   [volumes, nonEmptyTimeIndices, minTimePoints] = constructVolumes(smoothedHUImages);


    [numSlices, numTimePoints] = size(smoothedHUImages);
    nonEmptyTimeIndices = cell(numSlices, 1);
    
    for s = 1:numSlices
        nonEmptyTimeIndices{s} = find(~cellfun(@isempty, smoothedHUImages(s, :)));
    end
    
    minTimePoints = min(cellfun(@length, nonEmptyTimeIndices));
    if minTimePoints == 0
        error('At least one slice has no non-empty images.');
    end
    
    fprintf('Number of volumes to be constructed: %d\n', minTimePoints);
    volumes = cell(1, minTimePoints);
    
    for v = 1:minTimePoints
        [rows, cols] = size(smoothedHUImages{1, nonEmptyTimeIndices{1}(v)});
        volume_t = zeros(rows, cols, numSlices);
        for s = 1:numSlices
            currentImage = smoothedHUImages{s, nonEmptyTimeIndices{s}(v)};
            volume_t(:, :, s) = currentImage;
        end
        volumes{v} = volume_t;
    end
end
