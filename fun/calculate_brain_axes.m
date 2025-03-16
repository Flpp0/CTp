function calculate_brain_axes(baselineVolume, brainMask, outputDir)
% calculate_brain_axes Calculates the principal anatomical axes of the brain using PCA.
%
% Syntax:
%   calculate_brain_axes(baselineVolume, brainMask, outputDir)
%
% Description:
%   This function computes the principal axes of the brain starting from the
%   brain mask of the pre-contrast volume calculated by 'performSkullStripping' using Principal Component Analysis (PCA).
%   The computed principal components are classified into anatomical axes (Left-Right, Front-Back,
%   Top-Bottom). The function then computes the center of mass of the brain, saves the axes, PCA coefficients,
%   eigenvalues, and center of mass to a .mat file in the specified output directory,
%   and generates visualizations (with overlays) for each slice in the form of montages.
%
% Inputs:
%   baselineVolume - 3D array of the baseline CT image (e.g., [512 x 512 x numSlices]).
%   brainMask      - 3D binary mask of the brain (same dimensions as baselineVolume; 1 indicates brain, 0 background).
%   outputDir      - Directory path where the computed axes and visualizations will be saved.
%
% Example:
%   baselineVolume = load('baselineVolume.mat'); % Should load a 3D array
%   brainMask = load('brainMask.mat');           % Should load a binary 3D mask
%   outputDir = '/path/to/output';
%   calculate_brain_axes(baselineVolume, brainMask, outputDir);
%
% See also: classify_axes

    % Check
    if ~isequal(size(baselineVolume), size(brainMask))
        error('Dimensions of baselineVolume and brainMask must match.');
    end
    disp('Input validation passed. Dimensions are consistent.');

    % Extract Brain Mask Voxel Coordinates 
    [rows, cols, slices] = ind2sub(size(brainMask), find(brainMask));
    voxelCoords = [cols, rows, slices]; % [x (Left-Right), y (Front-Back), z (Top-Bottom)]
    disp(['Number of voxels in the brain mask: ', num2str(size(voxelCoords, 1))]);

    % Perform PCA 
    disp('Performing PCA on brain mask voxel coordinates...');
    [coeff, ~, latent] = pca(voxelCoords);
    disp('PCA completed.');

    % Classify Principal Components to Anatomical Axes
    axes = classify_axes(coeff);
    disp('Principal components classified to anatomical axes.');

    % Ensure Consistent Directionality Based on Anatomical Conventions
    % Left-Right: Left is negative X, so ensure X-axis points to the Right
    if axes.LeftRight(1) < 0
        axes.LeftRight = -axes.LeftRight;
    end

    % Front-Back: Front is positive Y, ensure Y-axis points to the Front
    if axes.FrontBack(2) < 0
        axes.FrontBack = -axes.FrontBack;
    end

    % Top-Bottom: Top is positive Z, ensure Z-axis points to the Top
    if axes.TopBottom(3) < 0
        axes.TopBottom = -axes.TopBottom;
    end

    centerOfMass = mean(voxelCoords, 1);

    % Save
    axesFile = fullfile(outputDir, 'principal_axes_corrected.mat');
    save(axesFile, 'axes', 'centerOfMass', 'coeff', 'latent');
    disp(['Principal axes and center of mass saved to: ', axesFile]);

    % Visualization
    disp('Visualizing principal axes using montages...');
    titles = {'Left-Right (LR) Plane', 'Front-Back (FB) Plane', 'Top-Bottom (TB) Plane'};
    colors = {'r', 'g', 'b'};  % Colors for each axis plane
    planes = {axes.LeftRight, axes.FrontBack, axes.TopBottom};
    numSlices = size(baselineVolume, 3); 

    for i = 1:3
        planeNormal = planes{i};
        d = dot(planeNormal, centerOfMass); % Plane constant
        
        % Preallocate 
        testFig = figure('Visible','off');
        colormap('gray');
        imagesc(baselineVolume(:, :, 1)); axis image;
        drawnow;
        testFrame = getframe(gca);
        [imH, imW, imC] = size(testFrame.cdata);
        close(testFig);
        axisImages = zeros(imH, imW, imC, numSlices, 'uint8');
        
        for sliceIdx = 1:numSlices
            fig = figure('Visible','off');
            colormap('gray');
            imagesc(baselineVolume(:, :, sliceIdx)); axis image; hold on;
            [X, Y] = meshgrid(1:size(baselineVolume,2), 1:size(baselineVolume,1));
            Z = sliceIdx * ones(size(X));
            planeMask = abs(planeNormal(1)*X + planeNormal(2)*Y + planeNormal(3)*Z - d) < 1;
            contour(planeMask, [0.5 0.5], colors{i}, 'LineWidth', 1.5);
            title(sprintf('Slice %d', sliceIdx));
            hold off;
            drawnow;
            frame = getframe(gca);
            axisImages(:,:,:,sliceIdx) = frame.cdata;
            close(fig);
        end

        % Montage
        figMontage = figure('Visible','off');
        montage(axisImages);
        title(titles{i}, 'FontSize', 16);
        
        % Save 
        montageFile = fullfile(outputDir, sprintf('%s_Montage.png', titles{i}));
        exportgraphics(gca, montageFile, 'Resolution', 300);
        close(figMontage);
        disp(['Saved montage for ', titles{i}, ' to: ', montageFile]);
    end

    disp('Principal axes visualization completed.');
end