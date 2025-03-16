function [common_time_base, allSlicesData, lastValidTimeSeconds, updatedBrainMasks] = calculateTCCs(registeredHUImages, smoothedHUImages, brainMasks, meanSlices, uniqueTimePointsSeconds, baseDir, patientCode, upsample_factor, axes, centerOfMass)
% calculateTCCs Calculates Tissue Concentration Curves (TCCs) with outlier removal.
%
% Syntax:
%   [common_time_base, allSlicesData, lastValidTimeSeconds, updatedBrainMasks] = calculateTCCs(registeredHUImages, smoothedHUImages, brainMasks, meanSlices, uniqueTimePointsSeconds, baseDir, patientCode, upsample_factor, axes, centerOfMass)
%
% Description:
%   This function computes a common time base and extracts TCCs for each slice from the registered
%   and smoothed images. It performs quadrant classification using the provided principal axes and
%   center of mass, then checks for outliers and updates brain masks.
%
% Inputs:
%   registeredHUImages - Cell array of registered images (nSlices x nTimePoints).
%   smoothedHUImages   - Cell array of smoothed images (nSlices x nTimePoints).
%   brainMasks         - 3D binary mask (nRows x nCols x nSlices) from skull stripping.
%   meanSlices         - 3D array of mean precontrast images.
%   uniqueTimePointsSeconds - Numeric vector of unique time points (in seconds).
%   baseDir            - Base directory for patient data.
%   patientCode        - Patient code used in folder names.
%   upsample_factor    - Upsampling factor.
%   axes               - Struct containing principal axes.
%   centerOfMass       - Center of mass coordinates.
%
% Outputs:
%   common_time_base   - Numeric vector of common time points for TCC interpolation.
%   allSlicesData      - Cell array where each element is a structure containing TCC data for a slice.
%   lastValidTimeSeconds - The maximum time point (in seconds) from the input.
%   updatedBrainMasks  - Updated 3D brain mask after outlier removal.
%
% Example:
%   [common_time_base, allSlicesData, lastValidTimeSeconds, updatedBrainMasks] = ...
%       calculateTCCs(registeredHUImages, smoothedHUImages, brainMasks, meanSlices, uniqueTimePointsSeconds, baseDir, patientCode, upsample_factor, axes, centerOfMass);
%
% See also: detectAndRemoveOutliers

    % Avoid using all the cores
    if isempty(gcp('nocreate'))
        numCores = feature('numcores');
        parpool('local', max(numCores - 1, 1));
    end

    % Saving paths
    tccDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)], 'TissueConcentrationCurves_Slices');
    outlierPlotsDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)], 'TCCsOutlierDetectionLinRegSlopeR2');
    
    if ~exist(tccDir, 'dir')
        mkdir(tccDir);
        disp(['Created TCC directory: ', tccDir]);
    end
    if ~exist(outlierPlotsDir, 'dir')
        mkdir(outlierPlotsDir);
        disp(['Created Outlier Plots directory: ', outlierPlotsDir]);
    end
    
    headerFile = fullfile(tccDir, 'TCC_header.mat');
    allSlicesFile = fullfile(tccDir, 'allSlicesData.mat');
    
    if exist(headerFile, 'file') && exist(allSlicesFile, 'file')
        disp('TCC data already exists. Loading precomputed data...');
        load(headerFile, 'common_time_base');
        load(allSlicesFile, 'allSlicesData');
        lastValidTimeSeconds = max(uniqueTimePointsSeconds);
    else
        disp('TCC data not found. Computing TCCs with outlier removal...');
        slopeThreshold = 0.4;
        r2Threshold = 0.6;
        plotOutliers = true;
        downsampleFactor = 10;
        uniqueTimePointsSecondsNorm = uniqueTimePointsSeconds - uniqueTimePointsSeconds(1);
        lastValidTimeSecondsNorm = max(uniqueTimePointsSecondsNorm);
        common_time_base = uniqueTimePointsSecondsNorm(1):1:lastValidTimeSecondsNorm;
        
        numSlices = size(registeredHUImages, 1);
        allSlicesData = cell(numSlices, 1);
        
        for sIdx = 1:numSlices
            fprintf('Processing slice: %d/%d\n', sIdx, numSlices);
            
            if ~any(brainMasks(:,:,sIdx))
                disp(['Skipping slice ', num2str(sIdx), ' due to empty brain mask.']);
                allSlicesData{sIdx} = [];
                continue;
            end
            
            sliceTimesIdx = find(~cellfun('isempty', registeredHUImages(sIdx, :)));
            if isempty(sliceTimesIdx)
                disp(['Skipping slice ', num2str(sIdx), ' due to empty time points.']);
                allSlicesData{sIdx} = [];
                continue;
            end
            
            sliceTimes = uniqueTimePointsSecondsNorm(sliceTimesIdx);
            [rows, cols] = size(registeredHUImages{sIdx, sliceTimesIdx(1)});
            voxIdx = find(brainMasks(:,:,sIdx));
            baselineSlice = meanSlices(:,:,sIdx);
            
            TCCs = zeros(length(voxIdx), length(common_time_base));
            for v = 1:length(voxIdx)
                [r, c] = ind2sub([rows, cols], voxIdx(v));
                intensity = zeros(1, length(sliceTimes));
                for tt = 1:length(sliceTimesIdx)
                    intensity(tt) = smoothedHUImages{sIdx, sliceTimesIdx(tt)}(r, c);
                end
                intensity = intensity - baselineSlice(r, c);
                TCCs(v, :) = interp1(sliceTimes, intensity, common_time_base, 'linear', 'extrap');
            end
            
            TCCs_original = TCCs;
            TCCs_filtered = movmean(TCCs_original, 3, 2);
            TCCs_filtered(TCCs_filtered < 0) = 0;
            
            % Quadrant Classification
            [voxelRows, voxelCols] = ind2sub([rows, cols], voxIdx);
            voxelPositions = [voxelCols, voxelRows, sIdx * ones(length(voxIdx), 1)];
            relativePositions = voxelPositions - centerOfMass;
            dotLR = relativePositions * axes.LeftRight;
            dotFB = relativePositions * axes.FrontBack;
            quadrant = zeros(length(voxIdx), 1);
            quadrant(dotLR > 0 & dotFB > 0) = 1;
            quadrant(dotLR < 0 & dotFB > 0) = 2;
            quadrant(dotLR < 0 & dotFB < 0) = 3;
            quadrant(dotLR > 0 & dotFB < 0) = 4;
            undefinedQuadrants = (quadrant == 0);
            if any(undefinedQuadrants)
                warning('Slice %d has voxels with undefined quadrants. Assigning default quadrant 1.', sIdx);
                quadrant(undefinedQuadrants) = 1;
            end
            
            sliceStruct = struct;
            sliceStruct.TCCs_original = TCCs_original;
            sliceStruct.TCCs_filtered = TCCs_filtered;
            sliceStruct.mask = brainMasks(:,:,sIdx);
            sliceStruct.voxIdx = voxIdx;
            sliceStruct.quadrant = quadrant;
            sliceStruct.rows = rows;
            sliceStruct.cols = cols;
            
            [sliceStruct, brainMasks(:,:,sIdx)] = detectAndRemoveOutliers(sliceStruct, brainMasks(:,:,sIdx), common_time_base, slopeThreshold, r2Threshold, outlierPlotsDir, sIdx, plotOutliers, downsampleFactor, baselineSlice);
            
            allSlicesData{sIdx} = sliceStruct;
            
            sliceFile = fullfile(tccDir, sprintf('slice_%d.mat', sIdx));
            save(sliceFile, 'sliceStruct', '-v7.3');
            disp(['Saved updated slice data for slice ', num2str(sIdx), ' to ', sliceFile]);
        end
        
        save(headerFile, 'common_time_base', '-v7.3');
        save(allSlicesFile, 'allSlicesData', '-v7.3');
        disp('TCC computation completed. Data saved successfully.');
        lastValidTimeSeconds = max(uniqueTimePointsSeconds);
    end
    
    updatedBrainMasks = brainMasks;

    %% Optional: Plot and Save TCC Subplots
    plotTCCs = true;
    if plotTCCs
        TCCPlotsDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)], 'TCCPlots');
        if ~exist(TCCPlotsDir, 'dir')
            mkdir(TCCPlotsDir);
            disp(['Created TCCPlots directory: ', TCCPlotsDir]);
        end
        
        maxSubplots = 36;
        downsampleFactor = 350; % Downsample factor for plotting TCCs
        numRowsPlot = 6;
        numColsPlot = 6;
        yMin = -20;
        yMax = 150;
        
        numSlices = length(allSlicesData);
        sliceFolders = cell(numSlices, 1);
        for sIdx = 1:numSlices
            sliceFolderName = sprintf('Slice_%d', sIdx);
            sliceFolders{sIdx} = fullfile(TCCPlotsDir, sliceFolderName);
            if ~exist(sliceFolders{sIdx}, 'dir')
                mkdir(sliceFolders{sIdx});
            end
        end
        
        parfor sIdx = 1:numSlices % Each slice is independent
            if isempty(allSlicesData{sIdx}) || ~isfield(allSlicesData{sIdx}, 'voxIdx') || isempty(allSlicesData{sIdx}.voxIdx)
                continue;
            end
            sliceData = allSlicesData{sIdx};
            TCCs_original = sliceData.TCCs_original;
            TCCs_filtered = sliceData.TCCs_filtered;
            voxIdx = sliceData.voxIdx;
            tccIndices = 1:downsampleFactor:length(voxIdx); % Downsample TCC indices
            nTCCs = length(tccIndices);
            if nTCCs == 0
                continue;
            end
            nFigures = ceil(nTCCs / maxSubplots);
            for figNum = 1:nFigures
                startIdx = (figNum - 1) * maxSubplots + 1;
                endIdx = min(figNum * maxSubplots, nTCCs);
                figName = sprintf('Slice_%d_Part_%d', sIdx, figNum);
                fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.8, 0.8]);
                t = tiledlayout(numRowsPlot, numColsPlot, 'TileSpacing', 'compact', 'Padding', 'compact');
                sgtitle(t, sprintf('Slice %d - TCCs (Downsampled by %d). Blue=Orig; Red=Filt.', sIdx, downsampleFactor), 'FontSize', 12, 'FontWeight', 'bold');
                for p = startIdx:endIdx
                    subplotIdx = p - (figNum - 1) * maxSubplots;
                    nexttile(subplotIdx);
                    rowIdx = tccIndices(p);
                    currentTCC_original = TCCs_original(rowIdx, :);
                    currentTCC_filtered = TCCs_filtered(rowIdx, :);
                    plot(common_time_base, currentTCC_original, 'b', 'LineWidth', 1.5);
                    hold on;
                    plot(common_time_base, currentTCC_filtered, 'r', 'LineWidth', 1.5);
                    hold off;
                    axis([common_time_base(1), common_time_base(end), yMin, yMax]);
                    grid on;
                    currentVoxIdx = voxIdx(rowIdx);
                    title(sprintf('Vox %d', currentVoxIdx), 'FontSize', 8);
                    if subplotIdx == 1 || subplotIdx == maxSubplots
                        xlabel('Time (s)');
                        ylabel('Intensity');
                    end
                end
                outFile = fullfile(sliceFolders{sIdx}, [figName, '.pdf']);
                try
                    exportgraphics(fig, outFile, 'ContentType', 'image', 'Resolution', 150);
                catch
                    print(fig, outFile, '-dpdf', '-bestfit');
                end
                close(fig);
                disp(['Saved TCC plot for ', figName, ' to ', outFile]);
            end
        end
        disp('All TCC plots created and saved as high-resolution PDFs.');
    end
end
