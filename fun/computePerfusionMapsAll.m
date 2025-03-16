function computePerfusionMapsAll(baseDir, patientCode, brainMasks, meanSlices, upsample_factor)
% computePerfusionMapsAll Computes and saves perfusion maps for sSVD, cSVD, and oSVD.
%
% Syntax:
%   computePerfusionMapsAll(baseDir, patientCode, brainMasks, meanSlices, upsample_factor)
%
% Description:
%   This function creates 3D volumes and computes perfusion maps (CBF, CBV, MTT, TMAX, and vessel mask)
%   for three methods: sSVD, cSVD, and oSVD. It loads a colormap from a CSV file and applies fixed color limits:
%       CBF:  [0, 35] % Custom
%       CBV:  [0, 4.5] % Custom
%       MTT:  [0, 12]
%       TMAX: [0, 15]
%   For each slice, the images are saved with the specified color limits and a colorbar.
%
% Inputs:
%   baseDir - Base directory for patient data.
%   patientCode - Patient code used in folder names.
%   brainMasks - 3D binary mask (nRows x nCols x nSlices).
%   meanSlices - 3D volume (nRows x nCols x nSlices) of mean precontrast images.
%   upsample_factor - Upsampling factor (1 means no upsampling).
%
% Example:
%   computePerfusionMapsAll(baseDir, patientCode, brainMasks, meanSlices, upsample_factor);
%
% See also: computePerfusionMaps, select_colormap

    % Define color limits for perfusion maps
    clim_cbf = [0, 35];
    clim_cbv = [0, 4.5];
    clim_mtt = [0, 12];
    clim_tmax = [0, 15];

    % Load Colormap
    clt_pma = readtable('PMA_lut.csv'); % ASIST Japan
    CLT = select_colormap(clt_pma, 'Siemens_CT'); % Comparison with Syngo.Via
    
    % sSVD Perfusion Maps
    disp('Computing sSVD Perfusion Maps...');
    PerfusionMapsDir_sSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
        ['Upsample_', num2str(upsample_factor)], 'PerfusionMaps_sSVD');
    if ~exist(PerfusionMapsDir_sSVD, 'dir')
        mkdir(PerfusionMapsDir_sSVD);
    end
    perfusionFile_sSVD = fullfile(PerfusionMapsDir_sSVD, 'PerfusionVolumes_sSVD.mat');
    
    if exist(perfusionFile_sSVD, 'file')
        disp('sSVD perfusion volumes already exist. Loading...');
        load(perfusionFile_sSVD, 'CBF_vol_sSVD','CBV_vol_sSVD','MTT_vol_sSVD','TMAX_vol_sSVD','vesselMask_vol_sSVD');
    else
        disp('No existing sSVD perfusion volumes found. Computing now...');
        CBF_Dir_sSVD = fullfile(PerfusionMapsDir_sSVD, 'CBF_Maps');
        CBV_Dir_sSVD = fullfile(PerfusionMapsDir_sSVD, 'CBV_Maps');
        MTT_Dir_sSVD = fullfile(PerfusionMapsDir_sSVD, 'MTT_Maps');
        TMAX_Dir_sSVD = fullfile(PerfusionMapsDir_sSVD, 'Tmax_Maps');
        VESSEL_Dir_sSVD = fullfile(PerfusionMapsDir_sSVD, 'VesselMask');
        if ~exist(CBF_Dir_sSVD, 'dir'), mkdir(CBF_Dir_sSVD); end
        if ~exist(CBV_Dir_sSVD, 'dir'), mkdir(CBV_Dir_sSVD); end
        if ~exist(MTT_Dir_sSVD, 'dir'), mkdir(MTT_Dir_sSVD); end
        if ~exist(TMAX_Dir_sSVD, 'dir'), mkdir(TMAX_Dir_sSVD); end
        if ~exist(VESSEL_Dir_sSVD, 'dir'), mkdir(VESSEL_Dir_sSVD); end
        
        [nRows, nCols, numSlices] = size(meanSlices);
        CBF_vol_sSVD   = zeros(nRows, nCols, numSlices);
        CBV_vol_sSVD   = zeros(nRows, nCols, numSlices);
        MTT_vol_sSVD   = zeros(nRows, nCols, numSlices);
        TMAX_vol_sSVD  = zeros(nRows, nCols, numSlices);
        vesselMask_vol_sSVD = false(nRows, nCols, numSlices);
        % Updated residue folder path:
        ResidueFunctionsDir_sSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
            ['Upsample_', num2str(upsample_factor)], 'ResidueFunctions', 'sSVD');
        
        for sIdx = 1:numSlices
            residueFile = fullfile(ResidueFunctionsDir_sSVD, sprintf('ResidueFunctions_Slice_%d_sSVD.mat', sIdx));
            if ~exist(residueFile, 'file')
                fprintf('No sSVD residue function file for slice %d. Skipping.\n', sIdx);
                continue;
            end
            load(residueFile, 'R_all','TCC_all','voxIdx','AIF','common_time_base');
            currentMask = brainMasks(:,:,sIdx);
            if isempty(currentMask) || ~any(currentMask(:))
                fprintf('No brain mask for slice %d. Skipping.\n', sIdx);
                continue;
            end
            [CBF_map, CBV_map, MTT_map, TMAX_map] = computePerfusionMaps(R_all, TCC_all, common_time_base, currentMask, AIF);
            vesselMask = (CBV_map > 8);
            CBF_map(vesselMask) = 0;
            CBV_map(vesselMask) = 0;
            MTT_map(vesselMask) = 0;
            TMAX_map(vesselMask) = 0;
            CBF_vol_sSVD(:,:,sIdx) = CBF_map;
            CBV_vol_sSVD(:,:,sIdx) = CBV_map;
            MTT_vol_sSVD(:,:,sIdx) = MTT_map;
            TMAX_vol_sSVD(:,:,sIdx) = TMAX_map;
            vesselMask_vol_sSVD(:,:,sIdx) = vesselMask;
            fprintf('Perfusion maps for slice %d (sSVD) computed.\n', sIdx);
            
            % Save individual slice images
            fig_cbf = figure('Visible','off');
            imagesc(CBF_map); axis image off; colormap(CLT); colorbar; clim(clim_cbf);
            title(sprintf('sSVD CBF (Slice %d)', sIdx));
            outFileCBF = fullfile(CBF_Dir_sSVD, sprintf('CBF_Slice_%d.png', sIdx));
            exportgraphics(fig_cbf, outFileCBF, 'Resolution',300);
            close(fig_cbf);
            
            fig_cbv = figure('Visible','off');
            imagesc(CBV_map); axis image off; colormap(CLT); colorbar; clim(clim_cbv);
            title(sprintf('sSVD CBV (Slice %d)', sIdx));
            outFileCBV = fullfile(CBV_Dir_sSVD, sprintf('CBV_Slice_%d.png', sIdx));
            exportgraphics(fig_cbv, outFileCBV, 'Resolution',300);
            close(fig_cbv);
            
            fig_mtt = figure('Visible','off');
            imagesc(MTT_map); axis image off; colormap(CLT); colorbar; clim(clim_mtt);
            title(sprintf('sSVD MTT (Slice %d)', sIdx));
            outFileMTT = fullfile(MTT_Dir_sSVD, sprintf('MTT_Slice_%d.png', sIdx));
            exportgraphics(fig_mtt, outFileMTT, 'Resolution',300);
            close(fig_mtt);
            
            fig_tmax = figure('Visible','off');
            imagesc(TMAX_map); axis image off; colormap(CLT); colorbar; clim(clim_tmax);
            title(sprintf('sSVD TMAX (Slice %d)', sIdx));
            outFileTmax = fullfile(TMAX_Dir_sSVD, sprintf('TMAX_Slice_%d.png', sIdx));
            exportgraphics(fig_tmax, outFileTmax, 'Resolution',300);
            close(fig_tmax);
        end
        save(perfusionFile_sSVD, 'CBF_vol_sSVD','CBV_vol_sSVD','MTT_vol_sSVD','TMAX_vol_sSVD','vesselMask_vol_sSVD','-v7.3');
        disp('Saved sSVD perfusion volumes.');
    end
    fig1 = figure; sliceViewer(CBF_vol_sSVD); colormap(CLT); title('Final CBF volume (sSVD)'); colorbar; clim(clim_cbf)
    fig2 = figure; sliceViewer(CBV_vol_sSVD); colormap(CLT); title('Final CBV volume (sSVD)'); colorbar; clim(clim_cbv)
    fig3 = figure; sliceViewer(MTT_vol_sSVD); colormap(CLT); title('Final MTT volume (sSVD)'); colorbar; clim(clim_mtt)
    fig4 = figure; sliceViewer(TMAX_vol_sSVD); colormap(CLT); title('Final TMAX volume (sSVD)'); colorbar; clim(clim_tmax)
    fig5 = figure; sliceViewer(vesselMask_vol_sSVD); colormap(gray); title('Final Vessel Mask (sSVD)'); colorbar;
    disp('sSVD perfusion map visualization complete.');
    
    % cSVD Perfusion Maps
    disp('Computing cSVD Perfusion Maps...');
    PerfusionMapsDir_cSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
        ['Upsample_', num2str(upsample_factor)], 'PerfusionMaps_cSVD');
    if ~exist(PerfusionMapsDir_cSVD, 'dir'), mkdir(PerfusionMapsDir_cSVD); end
    perfusionFile_cSVD = fullfile(PerfusionMapsDir_cSVD, 'PerfusionVolumes_cSVD.mat');
    
    if exist(perfusionFile_cSVD, 'file')
        disp('cSVD perfusion volumes already exist. Loading...');
        load(perfusionFile_cSVD, 'CBF_vol_cSVD','CBV_vol_cSVD','MTT_vol_cSVD','TMAX_vol_cSVD','vesselMask_vol_cSVD');
    else
        disp('No existing cSVD perfusion volumes found. Computing now...');
        CBF_Dir_cSVD = fullfile(PerfusionMapsDir_cSVD, 'CBF_Maps');
        CBV_Dir_cSVD = fullfile(PerfusionMapsDir_cSVD, 'CBV_Maps');
        MTT_Dir_cSVD = fullfile(PerfusionMapsDir_cSVD, 'MTT_Maps');
        TMAX_Dir_cSVD = fullfile(PerfusionMapsDir_cSVD, 'Tmax_Maps');
        VESSEL_Dir_cSVD = fullfile(PerfusionMapsDir_cSVD, 'VesselMask');
        if ~exist(CBF_Dir_cSVD, 'dir'), mkdir(CBF_Dir_cSVD); end
        if ~exist(CBV_Dir_cSVD, 'dir'), mkdir(CBV_Dir_cSVD); end
        if ~exist(MTT_Dir_cSVD, 'dir'), mkdir(MTT_Dir_cSVD); end
        if ~exist(TMAX_Dir_cSVD, 'dir'), mkdir(TMAX_Dir_cSVD); end
        if ~exist(VESSEL_Dir_cSVD, 'dir'), mkdir(VESSEL_Dir_cSVD); end
        
        [nRows, nCols, numSlices] = size(meanSlices);
        CBF_vol_cSVD   = zeros(nRows, nCols, numSlices);
        CBV_vol_cSVD   = zeros(nRows, nCols, numSlices);
        MTT_vol_cSVD   = zeros(nRows, nCols, numSlices);
        TMAX_vol_cSVD  = zeros(nRows, nCols, numSlices);
        vesselMask_vol_cSVD = false(nRows, nCols, numSlices);
        % Update residue folder path for cSVD:
        ResidueFunctionsDir_cSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
            ['Upsample_', num2str(upsample_factor)], 'ResidueFunctions', 'cSVD');
        
        for sIdx = 1:numSlices
            residueFile = fullfile(ResidueFunctionsDir_cSVD, sprintf('ResidueFunctions_Slice_%d_cSVD.mat', sIdx));
            if ~exist(residueFile, 'file')
                fprintf('No cSVD residue function file for slice %d. Skipping.\n', sIdx);
                continue;
            end
            load(residueFile, 'R_all','TCC_all','voxIdx','AIF','common_time_base');
            currentMask = brainMasks(:,:,sIdx);
            if isempty(currentMask) || ~any(currentMask(:))
                fprintf('No brain mask for slice %d. Skipping.\n', sIdx);
                continue;
            end
            [CBF_map, CBV_map, MTT_map, TMAX_map] = computePerfusionMaps(R_all, TCC_all, common_time_base, currentMask, AIF);
            vesselMask = (CBV_map > 8);
            CBF_map(vesselMask) = 0;
            CBV_map(vesselMask) = 0;
            MTT_map(vesselMask) = 0;
            TMAX_map(vesselMask) = 0;
            CBF_vol_cSVD(:,:,sIdx) = CBF_map;
            CBV_vol_cSVD(:,:,sIdx) = CBV_map;
            MTT_vol_cSVD(:,:,sIdx) = MTT_map;
            TMAX_vol_cSVD(:,:,sIdx) = TMAX_map;
            vesselMask_vol_cSVD(:,:,sIdx) = vesselMask;
            fprintf('Perfusion maps for slice %d (cSVD) computed.\n', sIdx);
            
            fig_cbv = figure('Visible','off');
            imagesc(CBV_map); axis image off; colormap(CLT); colorbar; clim(clim_cbv);
            title(sprintf('cSVD CBV (Slice %d)', sIdx));
            outFileCBV = fullfile(CBV_Dir_cSVD, sprintf('CBV_Slice_%d.png', sIdx));
            exportgraphics(fig_cbv, outFileCBV, 'Resolution',300);
            close(fig_cbv);
            
            fig_cbf = figure('Visible','off');
            imagesc(CBF_map); axis image off; colormap(CLT); colorbar; clim(clim_cbf);
            title(sprintf('cSVD CBF (Slice %d)', sIdx));
            outFileCBF = fullfile(CBF_Dir_cSVD, sprintf('CBF_Slice_%d.png', sIdx));
            exportgraphics(fig_cbf, outFileCBF, 'Resolution',300);
            close(fig_cbf);
            
            fig_mtt = figure('Visible','off');
            imagesc(MTT_map); axis image off; colormap(CLT); colorbar; clim(clim_mtt);
            title(sprintf('cSVD MTT (Slice %d)', sIdx));
            outFileMTT = fullfile(MTT_Dir_cSVD, sprintf('MTT_Slice_%d.png', sIdx));
            exportgraphics(fig_mtt, outFileMTT, 'Resolution',300);
            close(fig_mtt);
            
            fig_tmax = figure('Visible','off');
            imagesc(TMAX_map); axis image off; colormap(CLT); colorbar; clim(clim_tmax);
            title(sprintf('cSVD TMAX (Slice %d)', sIdx));
            outFileTmax = fullfile(TMAX_Dir_cSVD, sprintf('TMAX_Slice_%d.png', sIdx));
            exportgraphics(fig_tmax, outFileTmax, 'Resolution',300);
            close(fig_tmax);
        end
        save(perfusionFile_cSVD, 'CBF_vol_cSVD','CBV_vol_cSVD','MTT_vol_cSVD','TMAX_vol_cSVD','vesselMask_vol_cSVD','-v7.3');
        disp('Saved cSVD perfusion volumes.');
    end
    fig6 = figure; sliceViewer(CBF_vol_cSVD); colormap(CLT); title('Final CBF volume (cSVD)'); colorbar; clim(clim_cbf)
    fig7 = figure; sliceViewer(CBV_vol_cSVD); colormap(CLT); title('Final CBV volume (cSVD)'); colorbar; clim(clim_cbv)
    fig8 = figure; sliceViewer(MTT_vol_cSVD); colormap(CLT); title('Final MTT volume (cSVD)'); colorbar; clim(clim_mtt)
    fig9 = figure; sliceViewer(TMAX_vol_cSVD); colormap(CLT); title('Final TMAX volume (cSVD)'); colorbar; clim(clim_tmax)
    fig10 = figure; sliceViewer(vesselMask_vol_cSVD); colormap(gray); title('Final Vessel Mask (cSVD)'); colorbar;
    disp('cSVD perfusion map visualization complete.');
    
    % oSVD Perfusion Maps
    disp('Computing oSVD Perfusion Maps...');
    PerfusionMapsDir_oSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
        ['Upsample_', num2str(upsample_factor)], 'PerfusionMaps_oSVD');
    if ~exist(PerfusionMapsDir_oSVD, 'dir'), mkdir(PerfusionMapsDir_oSVD); end
    perfusionFile_oSVD = fullfile(PerfusionMapsDir_oSVD, 'PerfusionVolumes_oSVD.mat');
    
    if exist(perfusionFile_oSVD, 'file')
        disp('oSVD perfusion volumes already exist. Loading...');
        load(perfusionFile_oSVD, 'CBF_vol_oSVD','CBV_vol_oSVD','MTT_vol_oSVD','TMAX_vol_oSVD','vesselMask_vol_oSVD');
    else
        disp('No existing oSVD perfusion volumes found. Computing now...');
        CBF_Dir_oSVD = fullfile(PerfusionMapsDir_oSVD, 'CBF_Maps');
        CBV_Dir_oSVD = fullfile(PerfusionMapsDir_oSVD, 'CBV_Maps');
        MTT_Dir_oSVD = fullfile(PerfusionMapsDir_oSVD, 'MTT_Maps');
        TMAX_Dir_oSVD = fullfile(PerfusionMapsDir_oSVD, 'Tmax_Maps');
        VESSEL_Dir_oSVD = fullfile(PerfusionMapsDir_oSVD, 'VesselMask');
        if ~exist(CBF_Dir_oSVD, 'dir'), mkdir(CBF_Dir_oSVD); end
        if ~exist(CBV_Dir_oSVD, 'dir'), mkdir(CBV_Dir_oSVD); end
        if ~exist(MTT_Dir_oSVD, 'dir'), mkdir(MTT_Dir_oSVD); end
        if ~exist(TMAX_Dir_oSVD, 'dir'), mkdir(TMAX_Dir_oSVD); end
        if ~exist(VESSEL_Dir_oSVD, 'dir'), mkdir(VESSEL_Dir_oSVD); end
        
        [nRows, nCols, numSlices] = size(meanSlices);
        CBF_vol_oSVD   = zeros(nRows, nCols, numSlices);
        CBV_vol_oSVD   = zeros(nRows, nCols, numSlices);
        MTT_vol_oSVD   = zeros(nRows, nCols, numSlices);
        TMAX_vol_oSVD  = zeros(nRows, nCols, numSlices);
        vesselMask_vol_oSVD = false(nRows, nCols, numSlices);
        % Update residue folder path for oSVD:
        ResidueFunctionsDir_oSVD = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
            ['Upsample_', num2str(upsample_factor)], 'ResidueFunctions', 'oSVD');
        
        for sIdx = 1:numSlices
            residueFile = fullfile(ResidueFunctionsDir_oSVD, sprintf('ResidueFunctions_Slice_%d_oSVD.mat', sIdx));
            if ~exist(residueFile, 'file')
                fprintf('No oSVD residue function file for slice %d. Skipping.\n', sIdx);
                continue;
            end
            load(residueFile, 'R_all','TCC_all','OI_all','voxIdx','AIF','common_time_base');
            currentMask = brainMasks(:,:,sIdx);
            if isempty(currentMask) || ~any(currentMask(:))
                fprintf('No brain mask for slice %d. Skipping.\n', sIdx);
                continue;
            end
            [CBF_map, CBV_map, MTT_map, TMAX_map] = computePerfusionMaps(R_all, TCC_all, common_time_base, currentMask, AIF);
            vesselMask = (CBV_map > 8);
            CBF_map(vesselMask) = 0;
            CBV_map(vesselMask) = 0;
            MTT_map(vesselMask) = 0;
            TMAX_map(vesselMask) = 0;
            CBF_vol_oSVD(:,:,sIdx) = CBF_map;
            CBV_vol_oSVD(:,:,sIdx) = CBV_map;
            MTT_vol_oSVD(:,:,sIdx) = MTT_map;
            TMAX_vol_oSVD(:,:,sIdx) = TMAX_map;
            vesselMask_vol_oSVD(:,:,sIdx) = vesselMask;
            fprintf('Perfusion maps for slice %d (oSVD) computed.\n', sIdx);
            
            fig_cbf = figure('Visible','off');
            imagesc(CBF_map); axis image off; colormap(CLT); colorbar; clim(clim_cbf);
            title(sprintf('oSVD CBF (Slice %d)', sIdx));
            outFileCBF = fullfile(CBF_Dir_oSVD, sprintf('CBF_Slice_%d.png', sIdx));
            exportgraphics(fig_cbf, outFileCBF, 'Resolution',300);
            close(fig_cbf);
            
            fig_cbv = figure('Visible','off');
            imagesc(CBV_map); axis image off; colormap(CLT); colorbar; clim(clim_cbv);
            title(sprintf('oSVD CBV (Slice %d)', sIdx));
            outFileCBV = fullfile(CBV_Dir_oSVD, sprintf('CBV_Slice_%d.png', sIdx));
            exportgraphics(fig_cbv, outFileCBV, 'Resolution',300);
            close(fig_cbv);
            
            fig_mtt = figure('Visible','off');
            imagesc(MTT_map); axis image off; colormap(CLT); colorbar; clim(clim_mtt);
            title(sprintf('oSVD MTT (Slice %d)', sIdx));
            outFileMTT = fullfile(MTT_Dir_oSVD, sprintf('MTT_Slice_%d.png', sIdx));
            exportgraphics(fig_mtt, outFileMTT, 'Resolution',300);
            close(fig_mtt);
            
            fig_tmax = figure('Visible','off');
            imagesc(TMAX_map); axis image off; colormap(CLT); colorbar; clim(clim_tmax);
            title(sprintf('oSVD TMAX (Slice %d)', sIdx));
            outFileTmax = fullfile(TMAX_Dir_oSVD, sprintf('TMAX_Slice_%d.png', sIdx));
            exportgraphics(fig_tmax, outFileTmax, 'Resolution',300);
            close(fig_tmax);
        end
        save(perfusionFile_oSVD, 'CBF_vol_oSVD','CBV_vol_oSVD','MTT_vol_oSVD','TMAX_vol_oSVD','vesselMask_vol_oSVD','-v7.3');
        disp('Saved oSVD perfusion volumes.');
    end
    fig11 = figure; sliceViewer(CBF_vol_oSVD); colormap(CLT); title('Final CBF volume (oSVD)'); colorbar; clim(clim_cbf)
    fig12 = figure; sliceViewer(CBV_vol_oSVD); colormap(CLT); title('Final CBV volume (oSVD)'); colorbar; clim(clim_cbv)
    fig13 = figure; sliceViewer(MTT_vol_oSVD); colormap(CLT); title('Final MTT volume (oSVD)'); colorbar; clim(clim_mtt)
    fig14 = figure; sliceViewer(TMAX_vol_oSVD); colormap(CLT); title('Final TMAX volume (oSVD)'); colorbar; clim(clim_tmax)
    fig15 = figure; sliceViewer(vesselMask_vol_oSVD); colormap(gray); title('Final Vessel Mask (oSVD)'); colorbar;
    disp('oSVD perfusion map visualization complete.');
end
