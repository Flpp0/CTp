function [axes, centerOfMass, coeff, latent] = computePrincipalAxes(meanSlices, brainMasks, baseDir, patientCode, upsample_factor)
% computePrincipalAxes Computes or loads principal axes and the center of mass.
%
% Syntax:
%   [axes, centerOfMass, coeff, latent] = computePrincipalAxes(meanSlices, brainMasks, baseDir, patientCode, upsample_factor)
%
% Description:
%   This function checks for precomputed principal axes and center of mass stored in a file.
%   If not found, it calls the function calculate_brain_axes to compute
%   these values.
%
% Inputs:
%   meanSlices - 3D array of mean slices (nRows x nCols x nSlices).
%   brainMasks - 3D binary mask (nRows x nCols x nSlices).
%   baseDir - Base directory for patient data.
%   patientCode - Patient code used in folder names.
%   upsample_factor - Upsampling factor.
%
% Outputs:
%   axes - Struct containing the principal axes (e.g., LeftRight, FrontBack, TopBottom).
%   centerOfMass - Vector with center of mass coordinates.
%   coeff - PCA coefficient matrix.
%   latent - Eigenvalues from PCA.
%
% Example:
%   [axes, centerOfMass, coeff, latent] = computePrincipalAxes(meanSlices, brainMasks, baseDir, patientCode, upsample_factor);
%
% See also: calculate_brain_axes

    upsampleDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)]);
    outputDir = fullfile(upsampleDir, 'PrincipalAxes');
    if ~exist(outputDir, 'dir'), mkdir(outputDir); disp(['Created output directory: ', outputDir]); end
    axesFile = fullfile(outputDir, 'principal_axes_corrected.mat');
    
    if exist(axesFile, 'file')
        disp('Principal axes already calculated. Loading data...');
        load(axesFile, 'axes', 'centerOfMass', 'coeff', 'latent');
    else
        disp('Principal axes not found. Calculating and saving data...');
        calculate_brain_axes(meanSlices, brainMasks, outputDir);
        load(axesFile, 'axes', 'centerOfMass', 'coeff', 'latent');
    end
end

