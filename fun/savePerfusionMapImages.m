function savePerfusionMapImages(baseDir, patientCode, upsample_factor, meanSlices)
% savePerfusionMapImages Loads computed perfusion volumes and saves individual 
% slice images and montage images for each deconvolution method and perfusion parameter.
%
% Syntax:
%   savePerfusionMapImages(baseDir, patientCode, upsample_factor, meanSlices)
%
% Description:
%   This function loads the computed perfusion volumes for each deconvolution method
%   (sSVD, cSVD, and oSVD) from their respective .mat files stored in the upsample folder.
%   For each perfusion parameter (CBF, CBV, MTT, TMAX), it saves one image per slice 
%   (under "SliceImages") with a fixed colormap and color limits. In addition, using MATLAB's 
%   montage function, it creates a single montage figure that displays all slices together 
%   (under "SliceMontages") and saves that montage.
%
% Inputs:
%   baseDir         - Base directory for patient data.
%   patientCode     - Patient code used in folder names.
%   upsample_factor - Upsampling factor (e.g., 1 means no upsampling).
%   meanSlices      - 3D volume (nRows x nCols x nSlices) of mean precontrast images.
%
% See also: select_colormap

    % Do not use all the cores, leave 1 free
    if isempty(gcp('nocreate'))
        numCores = feature('numcores');
        parpool('local', max(numCores - 1, 1));
    end

    % Combinations
    methods = {'sSVD', 'cSVD', 'oSVD'};
    params = {'CBF', 'CBV', 'MTT', 'TMAX'};
    
    % Color limits: same as in computePerfusionMapsAll
    limits.CBF = [0, 35];
    limits.CBV = [0, 4.5];
    limits.MTT = [0, 12];
    limits.TMAX = [0, 15];
    
    % Colormap
    clt_pma = readtable('PMA_lut.csv'); % ASIST Japan
    CLT = select_colormap(clt_pma, 'Siemens_CT');
    
    % Loop
    for m = 1:length(methods)
        method = methods{m};
        % Path
        perfusionDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ...
            ['Upsample_', num2str(upsample_factor)], ['PerfusionMaps_', method]);
        switch method
            case 'sSVD'
                perfusionFile = fullfile(perfusionDir, 'PerfusionVolumes_sSVD.mat');
            case 'cSVD'
                perfusionFile = fullfile(perfusionDir, 'PerfusionVolumes_cSVD.mat');
            case 'oSVD'
                perfusionFile = fullfile(perfusionDir, 'PerfusionVolumes_oSVD.mat');
        end
        
        if ~exist(perfusionFile, 'file')
            fprintf('Perfusion file for %s does not exist. Skipping method %s.\n', method, method);
            continue;
        end
        
        % Perfusion volume data
        data = load(perfusionFile);
        
        CBF_vol = data.(['CBF_vol_', method]);
        CBV_vol = data.(['CBV_vol_', method]);
        MTT_vol = data.(['MTT_vol_', method]);
        TMAX_vol = data.(['TMAX_vol_', method]);
        
        % Create output folders for saving images:
        % For individual slice images:
        sliceImagesDir = fullfile(perfusionDir, 'SliceImages');
        if ~exist(sliceImagesDir, 'dir')
            mkdir(sliceImagesDir);
        end
        
        % Montage 
        montageDir = fullfile(perfusionDir, 'SliceMontages');
        if ~exist(montageDir, 'dir')
            mkdir(montageDir);
        end
        
        % Loop
        for p = 1:length(params)
            param = params{p};
            % Path
            paramSliceDir = fullfile(sliceImagesDir, param);
            if ~exist(paramSliceDir, 'dir')
                mkdir(paramSliceDir);
            end
            
            switch param
                case 'CBF'
                    vol = CBF_vol;
                case 'CBV'
                    vol = CBV_vol;
                case 'MTT'
                    vol = MTT_vol;
                case 'TMAX'
                    vol = TMAX_vol;
            end
            
            % Number of slices
            [~, ~, numSlices] = size(vol);
            
            % Save
            parfor sIdx = 1:numSlices
                fig = figure('Visible','off');
                imagesc(vol(:,:,sIdx));
                axis image off;
                colormap(CLT);
                colorbar;
                clim(limits.(param));
                title(sprintf('%s %s (Slice %d)', method, param, sIdx));
                outFile = fullfile(paramSliceDir, sprintf('%s_%s_Slice_%d.png', method, param, sIdx));
                exportgraphics(fig, outFile, 'Resolution',300);
                close(fig);
            end
            
            figMontage = figure('Visible','off');
            montage(vol, 'DisplayRange', limits.(param));
            colormap(CLT);
            colorbar;
            title(sprintf('Montage: %s %s', method, param));
            montageFile = fullfile(montageDir, sprintf('%s_%s_Montage.png', method, param));
            exportgraphics(gca, montageFile, 'Resolution',300);
            close(figMontage);
        end
        fprintf('Saved images and montages for method %s.\n', method);
    end
end
