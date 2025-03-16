function [registeredHUImages, volumeTransforms] = registerVolumes(interpolatedVolumes, uniqueSliceLocations, sortedPixelSpacing, baseDir, patientCode, upsample_factor)
% registerVolumes Performs 3D rigid registration on interpolated volumes using mutual information.
%
% Syntax:
%   [registeredHUImages, volumeTransforms] = registerVolumes(interpolatedVolumes, uniqueSliceLocations, sortedPixelSpacing, baseDir, patientCode, upsample_factor)
%
% Description:
%   This function either loads previously registered volumes from path or performs 3D rigid registration
%   on the provided full-resolution volumes using 'registerVolumes3D_MI_downsample' function.
%   The transformations are stored and the registered volumes are deconstructed back into a cell array of images.
%
% Inputs:
%   interpolatedVolumes - Cell array of upsampled volumes.
%   uniqueSliceLocations- Array of unique slice locations.
%   sortedPixelSpacing  - Cell array of pixel spacing values per slice.
%   baseDir             - Base directory for patient data.
%   patientCode         - Patient code (used for folder names).
%   upsample_factor     - Upsampling factor (1 means no upsampling).
%
% Outputs:
%   registeredHUImages  - Cell array (nSlices x nTimePoints) of registered images.
%   volumeTransforms    - Cell array of rigidtform3d objects for each volume.
%
% Example:
%   [regImages, trans] = registerVolumes(interpolatedVolumes, uniqueSliceLocations, sortedPixelSpacing, '/data', 'Patient01', 1);
%
% See also: registerVolumes3D_MI_downsample, imregtform, imwarp, rigidtform3d

    minTimePoints = length(interpolatedVolumes);
    registeredVolumesDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)], 'Volumes Registered');
    
    if exist(registeredVolumesDir, 'dir')
        disp('Loading existing registered volumes from individual files...');
        registeredVolumes = cell(1, minTimePoints);
        for v = 1:minTimePoints
            volumeDir = fullfile(registeredVolumesDir, sprintf('volume_%d', v));
            if ~exist(volumeDir, 'dir')
                error('Registered volume directory not found: %s', volumeDir);
            end
            sliceFiles = dir(fullfile(volumeDir, 'slice_*.mat'));
            if isempty(sliceFiles)
                error('No registered slice files found in volume directory %s.', volumeDir);
            end
            sliceNumbers = zeros(length(sliceFiles), 1);
            for k = 1:length(sliceFiles)
                tokens = regexp(sliceFiles(k).name, 'slice_(\d+)\.mat', 'tokens');
                if ~isempty(tokens)
                    sliceNumbers(k) = str2double(tokens{1}{1});
                else
                    error('Unexpected file name format: %s', sliceFiles(k).name);
                end
            end
            [~, idx] = sort(sliceNumbers);
            sliceFiles = sliceFiles(idx);
            for s = 1:length(sliceFiles)
                data = load(fullfile(volumeDir, sliceFiles(s).name), 'registered_slice');
                if s == 1
                    [rows, cols] = size(data.registered_slice);
                    numSlicesRegistered = length(sliceFiles);
                    registered_volume = zeros(rows, cols, numSlicesRegistered);
                end
                registered_volume(:, :, s) = data.registered_slice;
            end
            registeredVolumes{v} = registered_volume;
            fprintf('Loaded registered volume %d\n', v);
        end
        volumeTransformsFile = fullfile(registeredVolumesDir, 'volumeTransforms.mat');
        if exist(volumeTransformsFile, 'file')
            load(volumeTransformsFile, 'volumeTransforms');
        else
            error('volumeTransforms.mat not found in %s', registeredVolumesDir);
        end
        
        disp('Deconstructing registered volumes back into slices...');
        newNumSlices = size(registeredVolumes{1}, 3);
        registeredHUImages = cell(newNumSlices, minTimePoints);
        for v = 1:minTimePoints
            registeredVolume = registeredVolumes{v};
            for s = 1:newNumSlices
                registeredHUImages{s, v} = registeredVolume(:, :, s);
            end
        end
        disp('Registered volumes deconstructed into slices.');
    else
        disp('Performing 3D rigid registration on interpolated volumes using Mutual Information...');
        zThickness = abs(uniqueSliceLocations(2) - uniqueSliceLocations(1));
        pixelExtentInWorld = [sortedPixelSpacing{1}(1), sortedPixelSpacing{1}(2), zThickness];
        [registeredVolumes, volumeTransforms, executionTimes] = registerVolumes3D_MI_downsample(interpolatedVolumes);
        disp('3D rigid registration completed.');
        if ~exist(registeredVolumesDir, 'dir')
            mkdir(registeredVolumesDir);
        end
        for v = 1:minTimePoints
            registeredVolume = registeredVolumes{v};
            volumeDir = fullfile(registeredVolumesDir, sprintf('volume_%d', v));
            if ~exist(volumeDir, 'dir')
                mkdir(volumeDir);
            end
            numSlicesRegistered = size(registeredVolume, 3);
            for s = 1:numSlicesRegistered
                registered_slice = registeredVolume(:, :, s);
                save(fullfile(volumeDir, sprintf('slice_%d.mat', s)), 'registered_slice', '-v7.3');
            end
            fprintf('Registered volume %d slices saved.\n', v);
        end
        volumeTransformsFile = fullfile(registeredVolumesDir, 'volumeTransforms.mat');
        save(volumeTransformsFile, 'volumeTransforms', '-v7.3');
        disp('Deconstructing registered volumes back into slices...');
        newNumSlices = size(registeredVolumes{1}, 3);
        registeredHUImages = cell(newNumSlices, minTimePoints);
        for v = 1:minTimePoints
            registeredVolume = registeredVolumes{v};
            for s = 1:newNumSlices
                registeredHUImages{s, v} = registeredVolume(:, :, s);
            end
        end
        disp('Registered volumes deconstructed into slices.');
    end
end
