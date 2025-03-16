function interpolatedVolumes = upsampleVolumes(volumes, upsample_factor, baseDir, patientCode)
% upsampleVolumes Upsamples volumes in the Z-direction via spline interpolation.
%
% Syntax:
%   interpolatedVolumes = upsampleVolumes(volumes, upsample_factor, baseDir, patientCode)
%
% Description:
%   This function upsamples each 3D volume in the cell array along the Z-axis using spline or trilinear interpolation.
%   If a directory containing previously computed upsampled volumes exists, it loads them instead.
%   Otherwise, it computes and saves the upsampled volumes.
%
% Inputs:
%   volumes         - Cell array of 3D volumes.
%   upsample_factor - Upsampling factor (1 means no upsampling).
%   baseDir         - Base directory for saving results.
%   patientCode     - Patient code (used for constructing folder paths).
%
% Output:
%   interpolatedVolumes - Cell array of upsampled volumes.
%
% Example:
%   interpolatedVolumes = upsampleVolumes(volumes, 2, '/data', 'Patient01');
%
% See also: interp3, meshgrid

    minTimePoints = length(volumes);
    upsampledVolumesDir = fullfile(baseDir, patientCode, 'UpsampledVolumes', ['Upsample_', num2str(upsample_factor)]);
    
    if exist(upsampledVolumesDir, 'dir') && (numel(dir(upsampledVolumesDir)) > 2)
        disp('Loading existing upsampled volumes from individual files...');
        interpolatedVolumes = cell(1, minTimePoints);
        for v = 1:minTimePoints
            volumeDir = fullfile(upsampledVolumesDir, sprintf('volume_%d', v));
            if ~exist(volumeDir, 'dir')
                error('Volume directory not found: %s', volumeDir);
            end
            sliceFiles = dir(fullfile(volumeDir, 'slice_*.mat'));
            if isempty(sliceFiles)
                error('No upsampled slice files found in volume directory %s.', volumeDir);
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
                data = load(fullfile(volumeDir, sliceFiles(s).name), 'interpolated_slice');
                if s == 1
                    [rows, cols] = size(data.interpolated_slice);
                    numSlicesInterpolated = length(sliceFiles);
                    interpolated_volume = zeros(rows, cols, numSlicesInterpolated);
                end
                interpolated_volume(:, :, s) = data.interpolated_slice;
            end
            interpolatedVolumes{v} = interpolated_volume;
            disp(['Loaded upsampled volume ', num2str(v)]);
        end
        disp('All upsampled volumes loaded.');
    else
        mkdir(upsampledVolumesDir);
        if upsample_factor ~= 1
            disp('Upsampling volumes in Z via spline interpolation...');
            original_x = 1:size(volumes{1}, 2);
            original_y = 1:size(volumes{1}, 1);
            original_z = 1:size(volumes{1}, 3);
            query_x = original_x;
            query_y = original_y;
            query_z = 1:1/upsample_factor:original_z(end);
            [Xq, Yq, Zq] = meshgrid(query_x, query_y, query_z);
            interpolation_method = 'spline';
            for v = 1:length(volumes)
                current_volume = volumes{v};
                interpolated_volume = interp3(original_x, original_y, original_z, current_volume, Xq, Yq, Zq, interpolation_method);
                volumeDir = fullfile(upsampledVolumesDir, sprintf('volume_%d', v));
                if ~exist(volumeDir, 'dir')
                    mkdir(volumeDir);
                end
                numSlicesInterpolated = size(interpolated_volume, 3);
                for s = 1:numSlicesInterpolated
                    interpolated_slice = interpolated_volume(:, :, s);
                    save(fullfile(volumeDir, sprintf('slice_%d.mat', s)), 'interpolated_slice', '-v7.3');
                end
                disp(['Upsampled volume ', num2str(v), ' slices saved.']);
            end
            
            % Reload upsampled volumes
            interpolatedVolumes = cell(1, length(volumes));
            for v = 1:length(volumes)
                volumeDir = fullfile(upsampledVolumesDir, sprintf('volume_%d', v));
                if ~exist(volumeDir, 'dir')
                    error('Volume directory not found: %s', volumeDir);
                end
                sliceFiles = dir(fullfile(volumeDir, 'slice_*.mat'));
                if isempty(sliceFiles)
                    error('No upsampled slice files found in volume directory %s.', volumeDir);
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
                    data = load(fullfile(volumeDir, sliceFiles(s).name), 'interpolated_slice');
                    if s == 1
                        [rows, cols] = size(data.interpolated_slice);
                        numSlicesInterpolated = length(sliceFiles);
                        interpolated_volume = zeros(rows, cols, numSlicesInterpolated);
                    end
                    interpolated_volume(:, :, s) = data.interpolated_slice;
                end
                interpolatedVolumes{v} = interpolated_volume;
                disp(['Loaded upsampled volume ', num2str(v)]);
            end
            disp('All upsampled volumes loaded into interpolatedVolumes.');
        else
            disp('Upsampling factor is 1. Skipping upsampling.');
            interpolatedVolumes = volumes;
        end
    end
end
