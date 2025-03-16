function computeResidueFunctions(allSlicesData, common_time_base, AIF, upsample_factor, baseDir, patientCode)
% computeResidueFunctions Computes and saves residue functions using sSVD, cSVD, and oSVD.
%
% Syntax:
%   computeResidueFunctions(allSlicesData, common_time_base, AIF, upsample_factor, baseDir, patientCode)
%
% Description:
%   This function computes flow-scaled residue functions for each slice using three methods:
%     - Simple Truncated SVD (sSVD)
%     - Block-circulant cSVD: global threshold
%     - Block-circulant oSVD: oscillation index
%
%   For each slice and each method, the function first checks if the corresponding data file
%   already exists in the designated folder. If the file exists, the function skips computation.
%   Otherwise, the residue function is computed and saved.
%
%   Optionally, the function also creates subplots of the residue functions for each slice,
%   for each method, and saves them as high-resolution PDFs.
%
% Inputs:
%   allSlicesData    - Cell array containing TCC data for each slice.
%   common_time_base - Numeric vector of common time points.
%   AIF              - Arterial Input Function used for computing residue functions.
%   upsample_factor  - Upsampling factor.
%   baseDir          - Base directory for patient data.
%   patientCode      - Patient code used in folder names.
%
% See also: computeResidueFunctions_sSVD, computeResidueFunctions_cSVD, computeResidueFunctions_oSVD

numSlices = numel(allSlicesData);

%% sSVD
disp('Processing sSVD residue functions...');
ResidueFunctionsDir_sSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
    ['Upsample_', num2str(upsample_factor)], 'ResidueFunctions', 'sSVD');
if ~exist(ResidueFunctionsDir_sSVD, 'dir')
    mkdir(ResidueFunctionsDir_sSVD);
end
globalLambda = 0.2;
for sIdx = 1:numSlices
    sliceResidueFile = fullfile(ResidueFunctionsDir_sSVD, sprintf('ResidueFunctions_Slice_%d_sSVD.mat', sIdx));
    if exist(sliceResidueFile, 'file')
        fprintf('sSVD file exists for slice %d. Skipping computation.\n', sIdx);
        continue;
    end
    fprintf('Computing sSVD for slice %d...\n', sIdx);
    sliceData = allSlicesData{sIdx};
    if isempty(sliceData) || ~isfield(sliceData, 'TCCs_filtered') || isempty(sliceData.TCCs_filtered)
        disp('No data for this slice. Skipping.');
        continue;
    end
    TCCs_filtered = sliceData.TCCs_filtered;
    if size(TCCs_filtered,2) ~= length(common_time_base)
        error('Mismatch: TCC length does not match common_time_base length.');
    end
    R_all = computeResidueFunctions_sSVD(AIF, TCCs_filtered, globalLambda, common_time_base(2)-common_time_base(1));
    TCC_all = TCCs_filtered;
    voxIdx = sliceData.voxIdx;
    save(sliceResidueFile, 'R_all', 'TCC_all', 'voxIdx', 'AIF', 'common_time_base', 'globalLambda', '-v7.3');
    fprintf('sSVD residue functions saved for slice %d.\n', sIdx);
end
disp('sSVD residue functions processed.');

%% cSVD
disp('Processing cSVD residue functions...');
ResidueFunctionsDir_cSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
    ['Upsample_', num2str(upsample_factor)], 'ResidueFunctions', 'cSVD');
if ~exist(ResidueFunctionsDir_cSVD, 'dir')
    mkdir(ResidueFunctionsDir_cSVD);
end
lambdaRel = 0.1; m = 2;
for sIdx = 1:numSlices
    sliceResidueFile = fullfile(ResidueFunctionsDir_cSVD, sprintf('ResidueFunctions_Slice_%d_cSVD.mat', sIdx));
    if exist(sliceResidueFile, 'file')
        fprintf('cSVD file exists for slice %d. Skipping computation.\n', sIdx);
        continue;
    end
    fprintf('Computing cSVD for slice %d...\n', sIdx);
    sliceData = allSlicesData{sIdx};
    if isempty(sliceData) || ~isfield(sliceData, 'TCCs_filtered') || isempty(sliceData.TCCs_filtered)
        disp('No data for this slice. Skipping.');
        continue;
    end
    TCCs_filtered = sliceData.TCCs_filtered;
    if size(TCCs_filtered,2) ~= length(common_time_base)
        error('Mismatch: TCC length does not match common_time_base length.');
    end
    R_all = computeResidueFunctions_cSVD(AIF, TCCs_filtered, lambdaRel, m, common_time_base(2)-common_time_base(1));
    TCC_all = TCCs_filtered;
    voxIdx = sliceData.voxIdx;
    save(sliceResidueFile, 'R_all', 'TCC_all', 'voxIdx', 'AIF', 'common_time_base', 'lambdaRel', '-v7.3');
    fprintf('cSVD residue functions saved for slice %d.\n', sIdx);
end
disp('cSVD residue functions processed.');

%% oSVD
disp('Processing oSVD residue functions...');
ResidueFunctionsDir_oSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
    ['Upsample_', num2str(upsample_factor)], 'ResidueFunctions', 'oSVD');
if ~exist(ResidueFunctionsDir_oSVD, 'dir')
    mkdir(ResidueFunctionsDir_oSVD);
end
OI_limit = 0.06; m = 2; alphaCandidates = 0.01:0.005:0.95;
for sIdx = 1:numSlices
    sliceResidueFile = fullfile(ResidueFunctionsDir_oSVD, sprintf('ResidueFunctions_Slice_%d_oSVD.mat', sIdx));
    if exist(sliceResidueFile, 'file')
        fprintf('oSVD file exists for slice %d. Skipping computation.\n', sIdx);
        continue;
    end
    fprintf('Computing oSVD for slice %d...\n', sIdx);
    sliceData = allSlicesData{sIdx};
    if isempty(sliceData) || ~isfield(sliceData, 'TCCs_filtered') || isempty(sliceData.TCCs_filtered)
        disp('No data for this slice. Skipping.');
        continue;
    end
    TCCs_filtered = sliceData.TCCs_filtered;
    if size(TCCs_filtered,2) ~= length(common_time_base)
        error('Mismatch: TCC length does not match common_time_base length.');
    end
    [R_all, OI_all] = computeResidueFunctions_oSVD(AIF, TCCs_filtered, common_time_base(2)-common_time_base(1), OI_limit, m, alphaCandidates);
    TCC_all = TCCs_filtered;
    voxIdx = sliceData.voxIdx;
    save(sliceResidueFile, 'R_all', 'TCC_all', 'OI_all', 'voxIdx', 'AIF', 'common_time_base', 'OI_limit', 'alphaCandidates', '-v7.3');
    fprintf('oSVD residue functions saved for slice %d.\n', sIdx);
end
disp('oSVD residue functions processed.');

%% Optional: Plot and Save Residue Function Subplots
plotResidueFunctionsFlag = true;
if plotResidueFunctionsFlag
    % Avoid using all the cores for the parellelization
    if isempty(gcp('nocreate'))
        numCores = feature('numcores');
        parpool('local', max(numCores - 1, 1));
    end

    % Configuration parameters for the plotting (similar to the TCCs, better if they are the same)
    maxSubplots = 36;
    downsampleFactorPlot = 350; 
    numRowsPlot = 6;
    numColsPlot = 6;
    yMin = -0.01;  % Limits in the axis to compare the different subplots
    yMax = 0.03;
    
    % Directory
    methods = {'sSVD', 'cSVD', 'oSVD'};
    residueDirs = {ResidueFunctionsDir_sSVD, ResidueFunctionsDir_cSVD, ResidueFunctionsDir_oSVD};
    plotDirs = cell(size(methods));
    for mIdx = 1:length(methods)
        plotDirs{mIdx} = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
            ['Upsample_', num2str(upsample_factor)], 'ResiduePlots', methods{mIdx});
        if ~exist(plotDirs{mIdx}, 'dir')
            mkdir(plotDirs{mIdx});
        end
    end
    
    % For each method, plot residue functions for each slice in parallel
    parfor sIdx = 1:numSlices
        for mIdx = 1:length(methods)
            method = methods{mIdx};
            residueDir = residueDirs{mIdx};
            plotDir = plotDirs{mIdx};
            plotResidueFunctions(method, residueDir, plotDir, sIdx, maxSubplots, downsampleFactorPlot, numRowsPlot, numColsPlot, yMin, yMax);
        end
    end
    disp('Residue function plots created and saved as high-resolution PDFs.');
end
end

%% Helper Function: plotResidueFunctions
function plotResidueFunctions(method, residueDir, plotDir, sIdx, maxSubplots, downsampleFactor, numRowsPlot, numColsPlot, yMin, yMax)
    residueFile = fullfile(residueDir, sprintf('ResidueFunctions_Slice_%d_%s.mat', sIdx, method));
    if ~exist(residueFile, 'file')
        fprintf('No %s residue function file for slice %d. Skipping residue-plot.\n', method, sIdx);
        return;
    end
    % Load R_all ([V x T]), voxIdx, and common_time_base for plotting
    data = load(residueFile, 'R_all', 'voxIdx', 'common_time_base');
    if ~isfield(data, 'R_all') || isempty(data.R_all)
        fprintf('Empty %s residue for slice %d. Skipping.\n', method, sIdx);
        return;
    end
    R_all = data.R_all;
    voxIdx = data.voxIdx;
    common_time_base = data.common_time_base;
    
    nVox = size(R_all, 1);
    if nVox < 1
        fprintf('No voxels to plot for slice %d.\n', sIdx);
        return;
    end
    residueIndices = 1:downsampleFactor:nVox;
    nResidues = length(residueIndices);
    if nResidues == 0
        fprintf('No residues to plot after downsampling for slice %d.\n', sIdx);
        return;
    end
    nFigures = ceil(nResidues / maxSubplots);
    sliceFolder = fullfile(plotDir, sprintf('Slice_%d', sIdx));
    if ~exist(sliceFolder, 'dir')
        mkdir(sliceFolder);
    end
    for figNum = 1:nFigures
        startIdx = (figNum - 1) * maxSubplots + 1;
        endIdx = min(figNum * maxSubplots, nResidues);
        figName = sprintf('%s_Residues_Slice_%d_Part_%d', method, sIdx, figNum);
        fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
        t = tiledlayout(numRowsPlot, numColsPlot, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(t, sprintf('%s Residue Functions - Slice %d (Downsample=%d)', method, sIdx, downsampleFactor), 'FontSize', 12, 'FontWeight', 'bold');
        for p = startIdx:endIdx
            subplotIdx = p - (figNum - 1) * maxSubplots;
            nexttile(subplotIdx);
            rowIdx = residueIndices(p);
            oneResidue = R_all(rowIdx, :);
            plot(common_time_base, oneResidue, 'b', 'LineWidth', 1.5);
            axis([common_time_base(1), common_time_base(end), yMin, yMax]);
            grid on;
            title(sprintf('Vox %d', voxIdx(rowIdx)), 'FontSize', 8);
        end
        outFile = fullfile(sliceFolder, [figName, '.pdf']);
        try
            exportgraphics(fig, outFile, 'ContentType', 'image', 'Resolution', 150);
        catch
            print(fig, outFile, '-dpdf', '-bestfit');
        end
        close(fig);
        fprintf('Saved %s residue plot: %s to %s\n', method, figName, outFile);
    end
end
