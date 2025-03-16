function [updatedSliceStruct, updatedBrainMask] = detectAndRemoveOutliers(sliceStruct, brainMask, common_time_base, slopeThreshold, r2Threshold, outlierPlotsDir, sIdx, plotOutliers, downsampleFactor, baselineSlice)
    % detectAndRemoveOutliers Detects and removes outlier voxels from TCCs and updates the slice structure and brain mask.
    %
    % Syntax:
    %   [updatedSliceStruct, updatedBrainMask] = detectAndRemoveOutliers(sliceStruct, brainMask, common_time_base, slopeThreshold, r2Threshold, outlierPlotsDir, sIdx, plotOutliers, downsampleFactor, baselineSlice)
    %
    % Description:
    %   This function processes the TCCs in a given slice to
    %   identify outlier voxels using a simple linear regression model. Simple idea behind: we shall
    %   not have a constant drift that is increasing. For each voxel, a linear model is fit to its filtered TCC.
    %   Voxels are flagged as outliers if their absolute slope exceeds the specified slopeThreshold and their R² value 
    %   exceeds the r2Threshold (this ensure that we kind of have a constant drift).
    %   The outlier voxels are then removed from the TCC data, and the corresponding
    %   locations are marked as false in the brain mask. Optionally, plots of the outlier
    %   TCCs and an overlay of the outlier positions on the baseline slice are generated and saved.
    %
    % Inputs:
    %   sliceStruct      - Structure containing TCC data for the slice. Must include fields:
    %                      'TCCs_original', 'TCCs_filtered', 'voxIdx', and 'quadrant'.
    %   brainMask        - 2D binary mask for the slice.
    %   common_time_base - Numeric vector of common time points for TCC interpolation.
    %   slopeThreshold   - Threshold for the absolute slope in the linear regression.
    %   r2Threshold      - Threshold for the R² value in the linear regression.
    %   outlierPlotsDir  - Directory path to save outlier diagnostic plots.
    %   sIdx             - Slice index (used for labeling outputs).
    %   plotOutliers     - Boolean flag indicating whether to generate and save outlier plots.
    %   downsampleFactor - Factor by which to downsample the number of outlier plots.
    %   baselineSlice    - Baseline image slice used for overlaying outlier positions.
    %
    % Outputs:
    %   updatedSliceStruct - Updated slice structure with outlier voxels removed from TCC data.
    %   updatedBrainMask   - Updated binary brain mask with outlier voxel indices set to false.

    % Initialize 
    updatedSliceStruct = sliceStruct;
    updatedBrainMask = brainMask;

    if isempty(sliceStruct) || isempty(sliceStruct.voxIdx) || isempty(sliceStruct.TCCs_filtered)
        disp(['Skipping slice ', num2str(sIdx), ' due to empty slice data.']);
        return;
    end

    % Extract TCC Data
    TCCs_original = sliceStruct.TCCs_original; % Original TCCs
    TCCs_filtered = sliceStruct.TCCs_filtered; % Filtered TCCs
    numVoxels = size(TCCs_filtered, 1);
    numTimePoints = size(TCCs_filtered, 2);

    if numTimePoints < 2
        warning('Slice %d has less than 2 time points. Skipping outlier detection.', sIdx);
        return;
    end

    % Linear Regression
    X = [ones(length(common_time_base), 1), common_time_base(:)];
    beta = X \ TCCs_filtered'; % Regression coefficients [2 x numVoxels]
    TCCs_fitted = X * beta; % Fitted TCCs [numTimePoints x numVoxels]
    residuals = TCCs_filtered' - TCCs_fitted; % Residuals [numTimePoints x numVoxels]

    % R² Calculation
    ss_res = sum(residuals.^2, 1);
    ss_tot = sum((TCCs_filtered' - mean(TCCs_filtered', 1)).^2, 1);
    r2 = 1 - ss_res ./ ss_tot; % R² values [1 x numVoxels]

    % Slope Calculation
    slopes = beta(2, :)'; % Slopes [numVoxels x 1]

    % Identify Outliers (BOTH slope and R² criteria must be satisfied)
    isOutlier = (abs(slopes) > slopeThreshold) & (r2(:) > r2Threshold);
    numOutliers = sum(isOutlier);
    disp(['Number of outliers in slice ', num2str(sIdx), ': ', num2str(numOutliers)]);

    % Update Slice Struct and Brain Mask 
    if numOutliers > 0
        % Remove outliers from TCCs (both filtered and original)
        TCCs_original = TCCs_original(~isOutlier, :);
        TCCs_filtered = TCCs_filtered(~isOutlier, :);
        quadrant = sliceStruct.quadrant(~isOutlier); % <--- Update Quadrant Information

        % Update the slice structure
        updatedSliceStruct.TCCs_original = TCCs_original;
        updatedSliceStruct.TCCs_filtered = TCCs_filtered;
        updatedSliceStruct.quadrant = quadrant; % <--- Store Updated Quadrant

        % Remove outlier voxel indices
        outlierVoxelIndices = sliceStruct.voxIdx(isOutlier); % Indices of outliers
        updatedSliceStruct.voxIdx = sliceStruct.voxIdx(~isOutlier);

        % Update the brain mask
        updatedBrainMask(outlierVoxelIndices) = false;

        % Add the updated brain mask to the slice structure
        updatedSliceStruct.mask = updatedBrainMask;

        % Plot Outliers
        if plotOutliers
            downsampledOutlierIndices = find(isOutlier);
            downsampledOutlierIndices = downsampledOutlierIndices(1:downsampleFactor:end);
            numDownsampledOutliers = length(downsampledOutlierIndices);
            nFigures = ceil(numDownsampledOutliers / 36);

            sliceOutliersFolder = fullfile(outlierPlotsDir, sprintf('Slice_%d', sIdx));
            if ~exist(sliceOutliersFolder, 'dir')
                mkdir(sliceOutliersFolder);
            end

            for figNum = 1:nFigures
                startIdx = (figNum-1)*36 + 1;
                endIdx = min(figNum*36, numDownsampledOutliers);

                fig = figure('Visible', 'on', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
                tiledlayout(6, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
                sgtitle(sprintf('Slice %d - Outliers Part %d', sIdx, figNum));

                for p = startIdx:endIdx
                    subplotIdx = p - (figNum-1)*36;
                    nexttile(subplotIdx);

                    outlierIdx = downsampledOutlierIndices(p);
                    currentTCC_original = sliceStruct.TCCs_original(outlierIdx, :);
                    currentTCC_filtered = sliceStruct.TCCs_filtered(outlierIdx, :);
                    fittedLine = polyval(polyfit(common_time_base, currentTCC_filtered, 1), common_time_base);

                    % Plot
                    plot(common_time_base, currentTCC_original, 'b', 'LineWidth', 1.5);
                    hold on;
                    plot(common_time_base, currentTCC_filtered, 'r', 'LineWidth', 1.5);
                    plot(common_time_base, fittedLine, 'g--', 'LineWidth', 1.5);
                    hold off;
                    title(sprintf('Vox %d | Slope=%.2f | R²=%.2f', outlierIdx, slopes(outlierIdx), r2(outlierIdx)));
                end

                saveFile = fullfile(sliceOutliersFolder, sprintf('Outliers_Part_%d.pdf', figNum));
                exportgraphics(fig, saveFile, 'ContentType', 'vector');
                close(fig);
            end
        end

        % Plot Outlier Locations on Baseline Slice
        [outlierRows, outlierCols] = ind2sub(size(brainMask), outlierVoxelIndices);

        figOverlay = figure('Visible', 'on', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
        imshow(baselineSlice, []);
        hold on;
        scatter(outlierCols, outlierRows, 20, 'r', 'filled', 'MarkerEdgeColor', 'k');
        title(sprintf('Slice %d: Outlier Positions', sIdx), 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Column Index');
        ylabel('Row Index');
        hold off;

        overlayFile = fullfile(outlierPlotsDir, sprintf('Slice_%d_Outlier_Positions.pdf', sIdx));
        exportgraphics(figOverlay, overlayFile, 'ContentType', 'vector', 'Resolution', 300);
        close(figOverlay);
        disp(['Saved outlier positions overlay for slice ', num2str(sIdx), ' to ', overlayFile]);
    end
end
