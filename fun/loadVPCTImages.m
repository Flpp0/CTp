function [vpctImages, timePoints, sliceLocations, vpctPixelSpacing] = loadVPCTImages(info, dicomFiles, vpctDescription, pixelSpacing)
% loadVPCTImages Loads VPCT images and their PixelSpacing values.
%
% Syntax:
%   [vpctImages, timePoints, sliceLocations, vpctPixelSpacing] = loadVPCTImages(info, dicomFiles, vpctDescription, pixelSpacing)
%
% Description:
%   This function iterates through the DICOM info structures and loads the images
%   corresponding to the VPCT series. It extracts acquisition time, slice location,
%   and the associated PixelSpacing.
%
% Inputs:
%   info            - Cell array of DICOM info structures.
%   dicomFiles      - Structure array containing DICOM file information.
%   vpctDescription - String specifying the VPCT series description.
%   pixelSpacing    - Cell array of pixel spacing values.
%
% Outputs:
%   vpctImages      - Cell array of processed images (converted to HU).
%   timePoints      - Cell array of acquisition times.
%   sliceLocations  - Numeric array of slice locations.
%   vpctPixelSpacing- Cell array of PixelSpacing values for VPCT images.
%
% Example:
%   [vpctImages, timePoints, sliceLocations, vpctPixelSpacing] = loadVPCTImages(info, dicomFiles, vpctDescription, pixelSpacing);
%
% See also: dicomread, dicominfo

    numFiles = numel(info);
    timePoints = cell(numFiles, 1);
    sliceLocations = nan(numFiles, 1);
    vpctImages = cell(1, numFiles);
    vpctPixelSpacing = cell(1, numFiles);
    
    for i = 1:numFiles
        try
            if isfield(info{i}, 'SeriesDescription') && strcmp(info{i}.SeriesDescription, vpctDescription)
                acquisitionTime = num2str(info{i}.AcquisitionTime);
                sliceLocation = info{i}.SliceLocation;
                timePoints{i} = acquisitionTime;
                sliceLocations(i) = sliceLocation;
                vpctPixelSpacing{i} = pixelSpacing{i};
                if isfield(info{i}, 'RescaleSlope') && isfield(info{i}, 'RescaleIntercept')
                    rescaleSlope = info{i}.RescaleSlope;
                    rescaleIntercept = info{i}.RescaleIntercept;
                    % Use the original file name from the info structure if available
                    rawImage = double(dicomread(fullfile(info{i}.Filename)));
                    huImage = rawImage * rescaleSlope + rescaleIntercept;
                    vpctImages{i} = huImage;
                else
                    error('DICOM file missing RescaleSlope or RescaleIntercept');
                end
            else
                vpctImages{i} = [];
                vpctPixelSpacing{i} = [];
            end
        catch ME
            warning('Failed to read DICOM file: %s\nError: %s', dicomFiles(i).name, ME.message);
            vpctImages{i} = [];
            vpctPixelSpacing{i} = [];
        end
    end
    
    validIdx = ~cellfun('isempty', vpctImages);
    vpctImages = vpctImages(validIdx);
    timePoints = timePoints(validIdx);
    sliceLocations = sliceLocations(validIdx);
    vpctPixelSpacing = vpctPixelSpacing(validIdx);
end
