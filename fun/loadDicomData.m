function [info, images, pixelSpacing, uniqueSeriesDescriptions, dicomFiles] = loadDicomData(folderPath)
% loadDicomData Loads DICOM files from a folder and extracts metadata.
%
% Syntax:
%   [info, images, pixelSpacing, uniqueSeriesDescriptions, dicomFiles] = loadDicomData(folderPath)
%
% Description:
%   This function searches for all DICOM files in the specified folder,
%   loads their metadata using dicominfo, and extracts pixel spacing (if available)
%   and series descriptions. 
%
% Inputs:
%   folderPath - String containing the path to the folder with DICOM files.
%
% Outputs:
%   info                    - Cell array of DICOM info structures.
%   images                  - Cell array for image data.
%   pixelSpacing            - Cell array of pixel spacing values.
%   uniqueSeriesDescriptions- Cell array of unique series descriptions.
%   dicomFiles              - Structure array listing the DICOM files found.
%
% Example:
%   [info, images, pixelSpacing, uniqueSeriesDescriptions, dicomFiles] = loadDicomData('C:\Data\DICOM');
%
% See also: dicominfo

    dicomFiles = dir(fullfile(folderPath, '*.dcm'));
    numFiles = numel(dicomFiles);
    info = cell(1, numFiles);
    images = cell(1, numFiles);
    pixelSpacing = cell(1, numFiles);
    uniqueSeriesDescriptions = {};
    
    for i = 1:numFiles
        try
            info{i} = dicominfo(fullfile(folderPath, dicomFiles(i).name));
            if isfield(info{i}, 'SeriesDescription')
                currentDescription = info{i}.SeriesDescription;
                if ~any(strcmp(currentDescription, uniqueSeriesDescriptions))
                    uniqueSeriesDescriptions{end+1} = currentDescription;
                end
            end
            if isfield(info{i}, 'PixelSpacing')
                pixelSpacing{i} = info{i}.PixelSpacing;
            else
                warning('PixelSpacing not found for file: %s. Setting to [NaN, NaN].', dicomFiles(i).name);
                pixelSpacing{i} = [NaN, NaN];
            end
        catch ME
            warning('Failed to read DICOM info for file: %s\nError: %s', dicomFiles(i).name, ME.message);
            info{i} = [];
            pixelSpacing{i} = [NaN, NaN];
        end
    end
    
    validIdx = ~cellfun('isempty', info);
    info = info(validIdx);
    pixelSpacing = pixelSpacing(validIdx);
end

