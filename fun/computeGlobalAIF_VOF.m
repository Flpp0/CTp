function [aifStruct, vofStruct, globalAIFScaledVOF] = computeGlobalAIF_VOF(allSlicesData, common_time_base, meanSlices, brainMasks, axes, centerOfMass, baseDir, patientCode, upsample_factor)
% computeGlobalAIF_VOF Performs automatic global AIF and VOF detection and scales the AIF.
%
% Syntax:
%   [aifStruct, vofStruct, globalAIFScaledVOF] = computeGlobalAIF_VOF(allSlicesData, common_time_base, meanSlices, brainMasks, axes, centerOfMass, baseDir, patientCode, upsample_factor)
%
% Description:
%   This function calls computeGlobalAIF and computeGlobalVOF to determine the global arterial input function (AIF)
%   and venous output function (VOF), respectively, and then scales the AIF based on the VOF.
%
% Inputs:
%   allSlicesData   - Cell array containing TCC data for each slice.
%   common_time_base- Numeric vector of common time points.
%   meanSlices      - 3D array of mean slices.
%   brainMasks      - 3D binary mask.
%   axes            - Struct containing principal axes.
%   centerOfMass    - Center of mass coordinates.
%   baseDir         - Base directory.
%   patientCode     - Patient code used in folder names.
%   upsample_factor - Upsampling factor.
%
% Outputs:
%   aifStruct         - Struct containing the computed global AIF and related information.
%   vofStruct         - Struct containing the computed global VOF and related information.
%   globalAIFScaledVOF- The AIF scaled based on the VOF.
%
% See also: computeGlobalAIF, computeGlobalVOF

    saveDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)]);
    globalAIFDir = fullfile(saveDir, 'GlobalAIF');
    if ~exist(globalAIFDir, 'dir'), mkdir(globalAIFDir); disp(['Created directory for Global AIF: ', globalAIFDir]); end
    
    pAUC = 90; pIreg = 30; k = 5;
    aifStruct = computeGlobalAIF(allSlicesData, pAUC, pIreg, k, globalAIFDir, common_time_base, meanSlices, brainMasks);
    if isempty(aifStruct) || ~isfield(aifStruct, 'AIF')
        warning('No global AIF produced or missing AIF field.');
    else
        disp('Global AIF successfully computed.');
        disp(['Number of final AIF voxels: ', num2str(numel(aifStruct.voxelList))]);
        figure('Name','Final Global AIF'); plot(aifStruct.AIFTime, aifStruct.AIF, 'LineWidth',2);
        title('Final Global AIF (Unscaled)'); xlabel('Time (s)'); ylabel('Attenuation (HU)');
    end
    
    globalVOFDir = fullfile(saveDir, 'GlobalVOF');
    if ~exist(globalVOFDir, 'dir'), mkdir(globalVOFDir); disp(['Created directory for Global VOF: ', globalVOFDir]); end
    
    pAUC = 85; pIreg = 30; k = 5;
    vofStruct = computeGlobalVOF(allSlicesData, pAUC, pIreg, k, globalVOFDir, common_time_base, meanSlices, brainMasks, axes, centerOfMass);
    if isempty(vofStruct) || ~isfield(vofStruct, 'VOF')
        warning('No global VOF produced or missing VOF field.');
    else
        disp('Global VOF successfully computed.');
        disp(['Number of final VOF voxels: ', num2str(numel(vofStruct.voxelList))]);
        figure('Name','Final Global VOF','Visible','on'); plot(vofStruct.VOFTime, vofStruct.VOF, 'LineWidth',2);
        title('Final Global VOF (Unscaled)'); xlabel('Time (s)'); ylabel('Attenuation (HU)');
    end
    
    aucVOF = trapz(vofStruct.VOFTime, vofStruct.VOF);
    aucAIF = trapz(aifStruct.AIFTime, aifStruct.AIF);
    
    if aucVOF > 0
        globalAIFScaledVOF = aifStruct.AIF .* (aucVOF / aucAIF);
        disp('Plotting and saving scaled AIF, unscaled AIF, and VOF...');
        aifVofPlotsDir = fullfile(saveDir, 'AIF_VOF_Plots');
        if ~exist(aifVofPlotsDir, 'dir'), mkdir(aifVofPlotsDir); end
        [peakAIF, peakTimeAIFIdx] = max(aifStruct.AIF);
        peakTimeAIF = aifStruct.AIFTime(peakTimeAIFIdx);
        [peakVOF, peakTimeVOFIdx] = max(vofStruct.VOF);
        peakTimeVOF = vofStruct.VOFTime(peakTimeVOFIdx);
        [peakScaledAIF, peakTimeScaledAIFIdx] = max(globalAIFScaledVOF);
        peakTimeScaledAIF = aifStruct.AIFTime(peakTimeScaledAIFIdx);
        
        fig = figure('Name', 'AIF, Scaled AIF, and VOF', 'Position', [100, 100, 1000, 600]);
        hold on;
        % Unscaled AIF
        plot(aifStruct.AIFTime, aifStruct.AIF, 'Color', [1, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 3, 'DisplayName', 'Unscaled AIF');
        % Scaled AIF
        plot(aifStruct.AIFTime, globalAIFScaledVOF, 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 3, 'DisplayName', 'Scaled AIF (w.r.t VOF)');
        % VOF in blue
        plot(vofStruct.VOFTime, vofStruct.VOF, 'Color', [0, 0, 1], 'LineStyle', '-.', 'LineWidth', 3, 'DisplayName', 'VOF');
        % Peak markers
        plot(peakTimeAIF, peakAIF, 'o', 'MarkerEdgeColor', [1, 0.6, 0.6], 'MarkerFaceColor', [1, 0.6, 0.6], 'MarkerSize', 12, 'DisplayName', sprintf('AIF Peak (%.1f, %.1f)', peakTimeAIF, peakAIF));
        plot(peakTimeScaledAIF, peakScaledAIF, 'o', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'MarkerSize', 12, 'DisplayName', sprintf('Scaled AIF Peak (%.1f, %.1f)', peakTimeScaledAIF, peakScaledAIF));
        plot(peakTimeVOF, peakVOF, 'o', 'MarkerEdgeColor', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'MarkerSize', 12, 'DisplayName', sprintf('VOF Peak (%.1f, %.1f)', peakTimeVOF, peakVOF));
        xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Attenuation (HU)', 'FontSize', 16, 'FontWeight', 'bold');
        title('Unscaled AIF, Scaled AIF, and VOF with Peaks', 'FontSize', 18, 'FontWeight', 'bold');
        legend('show', 'Location', 'best', 'FontSize', 14);
        grid on; hold off;
        
        % PDF
        aifVofPlotFilePDF = fullfile(aifVofPlotsDir, 'AIF_VOF_ScaledAIF_Plot.pdf');
        exportgraphics(fig, aifVofPlotFilePDF, 'ContentType', 'vector', 'BackgroundColor', 'none');
        % PNG -> git hub readme
        aifVofPlotFilePNG = fullfile(aifVofPlotsDir, 'AIF_VOF_ScaledAIF_Plot.png');
        exportgraphics(fig, aifVofPlotFilePNG, 'Resolution',300);
        
        disp(['Plot saved as PDF: ', aifVofPlotFilePDF]);
        disp(['Plot saved as PNG: ', aifVofPlotFilePNG]);
    else
        warning('Integral of the VOF is negative; check AIF and VOF detection.');
        globalAIFScaledVOF = aifStruct.AIF; % fallback
    end
end
