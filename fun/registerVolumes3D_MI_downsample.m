function [registeredVolumes, volumeTransforms, executionTimes] = registerVolumes3D_MI_downsample(volumes, normalizationMethod, downsamplingFactor, maxIterationsFull)
% registerVolumes3D_MI_downsample Registers 3D volumes using a two-step 3D rigid registration with mutual information.
%
% Syntax:
%   [registeredVolumes, volumeTransforms, executionTimes] = registerVolumes3D_MI_downsample(volumes, normalizationMethod, downsamplingFactor, maxIterationsFull)
%
% Description:
%   This function performs 3D rigid registration on a series of volumes by employing a two-step process:
%     1. Downsample the XY plane of the input volumes for faster initial registration.
%     2. Use the transformations obtained from the downsampled data as initial estimates for full-resolution registration.
%   Mutual information is used as the similarity metric. The function also supports global intensity 
%   normalization if specified.
%
% Inputs:
%   volumes - Cell array of 3D volumes to be registered.
%   normalizationMethod - (Optional) String specifying the normalization method. Options are:
%                         'global' (to perform global normalization) or 
%                         'noNormalization' (to skip normalization). Default is 'nonormalization'.
%   downsamplingFactor - (Optional) Factor for downsampling on the XY plane (default is 2).
%   maxIterationsFull - (Optional) Maximum number of iterations for the full-resolution registration (default is 125).
%
% Outputs:
%   registeredVolumes - Cell array of full-resolution registered volumes.
%   volumeTransforms  - Cell array of rigidtform3d transformation objects corresponding to each volume.
%   executionTimes    - Structure containing execution times for the downsampled registration,
%                       full-resolution registration, and the total execution time.
%
% Example:
%   % Register volumes with default parameters:
%   [registeredVolumes, volumeTransforms, executionTimes] = registerVolumes3D_MI_downsample(volumes);
%
%   % Register volumes with global normalization, a downsampling factor of 2, and 150 full-resolution iterations:
%   [registeredVolumes, volumeTransforms, executionTimes] = registerVolumes3D_MI_downsample(volumes, 'global', 2, 150);
%
% See also: normalizeVolumeGlobal, downsampleXY

    % Overall time of execution. This function is the bottleneck of the whole code.
    % Parallization, when available, is used to speed up the calculations. 
    totalStartTime = tic; 

    % Default 
    if nargin < 2
        normalizationMethod = 'nonormalization';
    end
    if nargin < 3
        downsamplingFactor = 2;
    end
    if nargin < 4
        maxIterationsFull = 125;
    end

    % Convert volumes to single precision for performance 
    volumes = cellfun(@single, volumes, 'UniformOutput', false);

    % Normalization step (if requested)
    switch lower(normalizationMethod)
        case 'global'
            fprintf('Performing global normalization...\n');
            % Flatten volumes to compute global min and max
            allIntensities = cell2mat(volumes);
            globalMin = min(allIntensities(:));
            globalMax = max(allIntensities(:));
            % Normalize each volume
            volumes = cellfun(@(v) normalizeVolumeGlobal(v, globalMin, globalMax), volumes, 'UniformOutput', false);
        case 'nonormalization'
            fprintf('Skipping normalization (using raw intensity values)...\n');
        otherwise
            error('Unsupported normalization method. Use ''global'' or ''noNormalization''.');
    end

    % Step 1: Downsample the volumes in the XY plane
    fprintf('Downsampling volumes by a factor of %d on the XY plane...\n', downsamplingFactor);
    downsampledVolumes = cellfun(@(v) downsampleXY(v, downsamplingFactor), volumes, 'UniformOutput', false);

    % Initialize transformations for downsampled registration
    numVolumes = length(volumes);
    downsampledTransforms = cell(1, numVolumes);
    fixedVolumeDownsampled = downsampledVolumes{1};
    downsampledTransforms{1} = rigidtform3d(); % Identity transform for the first volume

    % Configure optimizer and metric for downsampled registration
    [optimizer, metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 500;
    optimizer.InitialRadius = 5e-4;
    optimizer.Epsilon = 1.5e-5;

    % Set up progress logging for downsampled registration
    progressLogDownsampled = parallel.pool.DataQueue;
    afterEach(progressLogDownsampled, @(msg) fprintf('%s\n', msg));

    % Time downsampled registration
    downsampledStartTime = tic;

    % Register downsampled volumes (parallelized)
    fprintf('Performing initial registration on downsampled volumes...\n');
    parfor t = 2:numVolumes
        try
            movingVolumeDownsampled = downsampledVolumes{t};
            tform = imregtform(movingVolumeDownsampled, fixedVolumeDownsampled, 'rigid', optimizer, metric);
            downsampledTransforms{t} = tform;
            send(progressLogDownsampled, sprintf('Downsampled volume %d registered successfully.', t));
        catch ME
            warning('Registration failed for downsampled volume %d. Error: %s', t, ME.message);
            downsampledTransforms{t} = rigidtform3d();
            send(progressLogDownsampled, sprintf('Registration failed for downsampled volume %d.', t));
        end
    end

    executionTimes.downsampledRegistration = toc(downsampledStartTime);

    % Step 2: Full-resolution registration refinement
    fprintf('Refining registration on full-resolution volumes...\n');
    tempRegisteredVolumes = cell(1, numVolumes);
    tempVolumeTransforms = cell(1, numVolumes);
    fixedVolume = volumes{1};
    tempRegisteredVolumes{1} = fixedVolume;
    tempVolumeTransforms{1} = rigidtform3d();

    % Reconfigure optimizer for full resolution
    % The following registration parameters can be fine tuned
    optimizer.MaximumIterations = maxIterationsFull;
    optimizer.InitialRadius = 2.5e-4;
    optimizer.Epsilon = 1e-6;

    progressLogFull = parallel.pool.DataQueue;
    afterEach(progressLogFull, @(msg) fprintf('%s\n', msg));
    fullResolutionStartTime = tic;

    parfor t = 2:numVolumes
        try
            movingVolume = volumes{t};
            initialTransform = downsampledTransforms{t};
            tform = imregtform(movingVolume, fixedVolume, 'rigid', optimizer, metric, 'InitialTransformation', initialTransform);
            registeredVolume = imwarp(movingVolume, tform, 'OutputView', imref3d(size(fixedVolume)), 'InterpolationMethod', 'linear');
            tempRegisteredVolumes{t} = registeredVolume;
            tempVolumeTransforms{t} = tform;
            send(progressLogFull, sprintf('Full-resolution volume %d registered successfully.', t));
        catch ME
            warning('Full-resolution registration failed for volume %d. Error: %s', t, ME.message);
            tempRegisteredVolumes{t} = volumes{t};
            tempVolumeTransforms{t} = rigidtform3d();
            send(progressLogFull, sprintf('Registration failed for full-resolution volume %d.', t));
        end
    end

    executionTimes.fullResolutionRegistration = toc(fullResolutionStartTime);
    registeredVolumes = tempRegisteredVolumes;
    volumeTransforms = tempVolumeTransforms;
    executionTimes.total = toc(totalStartTime);
end

