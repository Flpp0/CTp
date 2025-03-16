function aifStruct = computeGlobalAIF(allSlicesData, pAUC, pIreg, k, outDir, common_time_base, meanSlices, brainMasks)
% computeGlobalAIF Computes the global Arterial Input Function (AIF) from filtered tissue concentration curves.
%
% Syntax:
%   aifStruct = computeGlobalAIF(allSlicesData, pAUC, pIreg, k, outDir, common_time_base, meanSlices, brainMasks)
%
% Description:
%   This function aggregates filtered time-attenuation curves (TCCs) from multiple slices,
%   applies "second AUC" filtering (keeping only the top (100-pAUC)% of curves) and roughness
%   filtering (removing the top pIreg% of rough curves based on the second derivative norm),
%   then performs a two-step k-means clustering on the normalized curves. The final cluster (with
%   the lowest first moment) is mapped back to the original raw TACs amplitudes to yield the global AIF.
%
%   The automatic selection of the arterial input function is based on a clustering approach
%   described by Mouridsen et al. (2006) with adapted numbers:
%
%       K. Mouridsen, S. Christensen, L. Gyldensted, and L. Østergaard, “Automatic selection of arterial 
%       input function using cluster analysis,” Magnetic Resonance in Med, vol. 55, no. 3, pp. 524–531, 
%       Mar. 2006, doi: 10.1002/mrm.20759.
%
%   Diagnostic plots are generated and saved in the specified output directory.
%
% Inputs:
%   allSlicesData    - Cell array (one element per slice) with fields:
%                        .TCCs_filtered (V x T): Filtered time-attenuation curves.
%                        .voxIdx        (V x 1): Voxel indices in the slice.
%                        .rows, .cols   : Dimensions of the slice.
%                        .mask          : Binary mask (nRows x nCols).
%   pAUC             - Percentage threshold (0 to 100) for "second AUC" filtering.
%   pIreg            - Percentage threshold (0 to 100) for roughness filtering.
%   k                - Number of clusters for k-means clustering.
%   outDir           - Directory for saving diagnostic plots.
%   common_time_base - 1 x T vector of common time points.
%   meanSlices       - 3D array (nRows x nCols x numSlices) of mean precontrast images.
%   brainMasks       - 3D logical array (nRows x nCols x numSlices) representing brain masks.
%
% Output:
%   aifStruct - Structure containing:
%                 .AIFTime   - The common time base.
%                 .AIF       - The global AIF curve (raw TAC mean of the final cluster).
%                 .voxelList - Global indices (in the concatenated TCC array) of the final cluster.
%                 .Xcoord    - Row coordinates of the final cluster voxels.
%                 .Ycoord    - Column coordinates of the final cluster voxels.
%                 .sliceIdx  - Slice indices for each voxel.
%
% See also: computeRoughness, computeFirstMoment

    %% 1) Gather TCCs from all slices into a big array 
    C_original_all = [];
    Xcoord_all     = [];
    Ycoord_all     = [];
    sliceIdx_all   = [];
    bigIdx_all     = [];
    offset         = 0;
    numSlices = numel(allSlicesData);

    for sIdx = 1:numSlices
        sliceData = allSlicesData{sIdx};
        if isempty(sliceData) || ~isfield(sliceData, 'TCCs_filtered') || isempty(sliceData.TCCs_filtered)
            continue;
        end
        TACs   = sliceData.TCCs_filtered;
        voxIdx = sliceData.voxIdx;
        rows   = sliceData.rows;
        cols   = sliceData.cols;
        [Xc, Yc] = ind2sub([rows, cols], voxIdx);
        nV = size(TACs,1);
        C_original_all = [C_original_all; TACs];
        newGlobalIdx = (offset + 1 : offset + nV).';
        bigIdx_all   = [bigIdx_all; newGlobalIdx];
        Xcoord_all   = [Xcoord_all; Xc];
        Ycoord_all   = [Ycoord_all; Yc];
        sliceIdx_tmp = repmat(sIdx, [nV, 1]);
        sliceIdx_all = [sliceIdx_all; sliceIdx_tmp];
        offset = offset + nV;
    end

    if isempty(C_original_all)
        warning('No TCCs found across all slices. Returning empty struct.');
        aifStruct = [];
        return;
    end

    C_original  = C_original_all;
    Xcoord      = Xcoord_all;
    Ycoord      = Ycoord_all;
    sliceIdxV   = sliceIdx_all;
    bigIdx      = bigIdx_all;
    [T] = size(C_original,2);
    if T ~= length(common_time_base)
        warning('C_original time dimension mismatch with common_time_base.');
    end

    %% 2) "Second AUC" Filtering
    AUC = trapz(common_time_base, C_original, 2);
    Nvox = size(C_original,1);
    [~, idxAUCsort] = sort(AUC, 'ascend');
    cutAUC = round(Nvox * (pAUC / 100));
    if cutAUC < Nvox
        idxKeepAUC = idxAUCsort(cutAUC+1:end);
    else
        warning('No curves remain after second AUC filtering.');
        aifStruct = [];
        return;
    end
    C         = C_original(idxKeepAUC, :);
    Xcoord    = Xcoord(idxKeepAUC);
    Ycoord    = Ycoord(idxKeepAUC);
    sliceIdxV = sliceIdxV(idxKeepAUC);
    bigIdx    = bigIdx(idxKeepAUC);
    AUC_sel   = AUC(idxKeepAUC);

    if isempty(C)
        warning('No curves remain after second AUC filtering...');
        aifStruct = [];
        return;
    end

    % Optional AUC histogram
    figAUC = figure('Visible','off');
    histogram(AUC_sel, 100);
    title(sprintf('AUC distribution after pAUC=%.0f filtering', pAUC));
    xlabel('AUC'); ylabel('# Voxels');
    saveas(figAUC, fullfile(outDir, 'GlobalAIF_AUC_histogram.png'));
    close(figAUC);

    %% 3) Normalize by AUC for further filtering
    areaScale   = trapz(common_time_base, C, 2);
    Cunit_area  = C ./ areaScale;

    %% 4) Roughness Filtering 
    roughVals = computeRoughness(Cunit_area, common_time_base);
    [~, idxRsort] = sort(roughVals, 'descend');
    cutR = round(length(idxRsort)*(pIreg/100));
    if cutR < length(idxRsort)
        idxKeepR = idxRsort(cutR+1:end);
    else
        warning('No curves remain after roughness filtering.');
        aifStruct = [];
        return;
    end

    Cunit_area = Cunit_area(idxKeepR, :);
    C          = C(idxKeepR, :);
    Xcoord     = Xcoord(idxKeepR);
    Ycoord     = Ycoord(idxKeepR);
    sliceIdxV  = sliceIdxV(idxKeepR);
    bigIdx     = bigIdx(idxKeepR);

    if isempty(Cunit_area)
        warning('No curves remain after roughness filtering...');
        aifStruct = [];
        return;
    end

    % Optional roughness histogram
    figRough = figure('Visible','off');
    histogram(roughVals(idxKeepR), 100);
    title(sprintf('Roughness after removing top %.0f%%', pIreg));
    xlabel('Roughness'); ylabel('# Voxels');
    saveas(figRough, fullfile(outDir, 'GlobalAIF_Roughness_histogram.png'));
    close(figRough);

    %%  5) First Clustering (k-means)
    [clusterIdx1, Cmeans1] = kmeans(C, k, 'Replicates', 5, 'MaxIter', 1500);
    fm1 = computeFirstMoment(Cmeans1, common_time_base);
    [~, minIdx1] = min(fm1);

    cmap_full = autumn(256);
    figClust1 = figure('Visible','off');
    hold on;
    fm_min = min(fm1);
    fm_max = max(fm1);
    if fm_min == fm_max
        idxColor1 = repmat(128, [k,1]);
    else
        idxColor1 = round((fm1 - fm_min)/(fm_max - fm_min)*255)+1;
    end
    for kk = 1:k
        curve = Cmeans1(kk,:);
        colorNow = cmap_full(idxColor1(kk), :);
        if kk==minIdx1
            plot(common_time_base, curve, 'Color','k','LineWidth',3);
        else
            plot(common_time_base, curve, 'Color', colorNow, 'LineWidth',2);
        end
    end
    title('First Clustering (Global AIF) - Mean Curves');
    xlabel('Time (s)'); ylabel('Raw Attenuation');
    hold off;
    saveas(figClust1, fullfile(outDir, 'GlobalAIF_FirstClustering_MeanCurves.png'));
    close(figClust1);

    % Subset for minimal first-moment cluster
    idxSubset1 = find(clusterIdx1 == minIdx1);
    Csubset1       = Cunit_area(idxSubset1, :);
    CrawSubset1    = C(idxSubset1, :);
    bigIdx_sub1    = bigIdx(idxSubset1);
    Xsub1          = Xcoord(idxSubset1);
    Ysub1          = Ycoord(idxSubset1);
    Ssub1          = sliceIdxV(idxSubset1);

    %% 6) Second Clustering on that subset
    [clusterIdx2, Cmeans2] = kmeans(CrawSubset1, k, 'Replicates',5,'MaxIter',1000);
    fm2 = computeFirstMoment(Cmeans2, common_time_base);
    [~, minIdx2] = min(fm2);

    figClust2 = figure('Visible','off');
    hold on;
    fm2_min = min(fm2);
    fm2_max = max(fm2);
    if fm2_min == fm2_max
        idxColor2 = repmat(128, [k,1]);
    else
        idxColor2 = round((fm2 - fm2_min)/(fm2_max - fm2_min)*255)+1;
    end
    for kk = 1:k
        curve = Cmeans2(kk,:);
        colorC = cmap_full(idxColor2(kk), :);
        if kk==minIdx2
            plot(common_time_base, curve, 'Color','k','LineWidth',3);
        else
            plot(common_time_base, curve, 'Color', colorC, 'LineWidth',2);
        end
    end
    title('Second Clustering (Global AIF) - Mean Curves');
    xlabel('Time (s)'); ylabel('Raw Attenuation');
    hold off;
    saveas(figClust2, fullfile(outDir, 'GlobalAIF_SecondClustering_MeanCurves.png'));
    close(figClust2);

    % Identify final subset from second clustering
    idxFinal2       = find(clusterIdx2 == minIdx2);
    finalIdx_all    = bigIdx_sub1(idxFinal2);
    Xfinal          = Xsub1(idxFinal2);
    Yfinal          = Ysub1(idxFinal2);
    Sfinal          = Ssub1(idxFinal2);

    finalClusterOriginal = C_original_all(finalIdx_all, :);
    AIFraw_unscaled      = mean(finalClusterOriginal, 1);
    AIFtime              = common_time_base;

    % Build Output Structure
    aifStruct = struct();
    aifStruct.AIFTime    = AIFtime;
    aifStruct.AIF        = AIFraw_unscaled;   
    aifStruct.voxelList  = finalIdx_all;
    aifStruct.Xcoord     = Xfinal;
    aifStruct.Ycoord     = Yfinal;
    aifStruct.sliceIdx   = Sfinal;

    fprintf('computeGlobalAIF: final cluster has %d voxels.\n', numel(finalIdx_all));

    %% 7) Overlays per slice
    overlayDir = fullfile(outDir, 'GlobalAIF_VoxelOverlays');
    if ~exist(overlayDir, 'dir'), mkdir(overlayDir); end

    [nRows, nCols, nSlcs] = size(meanSlices);
    for s = 1:nSlcs
        thisSliceMask = (Sfinal == s);
        if ~any(thisSliceMask), continue; end
        Xsel = Xfinal(thisSliceMask);
        Ysel = Yfinal(thisSliceMask);
        figName = sprintf('Slice_%d_GlobalAIF_voxels', s);
        fig = figure('Visible','on');
        imshow(meanSlices(:,:,s), []);
        hold on;
        scatter(Ysel, Xsel, 20, 'r','filled');
        title(sprintf('Global AIF Voxels - Slice %d', s));
        outFile = fullfile(overlayDir, [figName, '.png']);
        saveas(fig, outFile);
        close(fig);
        fprintf('Saved overlay: %s\n', outFile);
    end

    % Save a montage of all slices with marked voxels 
    allOverlayImages = zeros(nRows, nCols, 3, nSlcs, 'uint8');
    for s = 1:nSlcs
        fig = figure('Visible','off');
        imshow(meanSlices(:,:,s), []);
        hold on;
        thisSliceMask = (Sfinal == s);
        if any(thisSliceMask)
            Xsel = Xfinal(thisSliceMask);
            Ysel = Yfinal(thisSliceMask);
            scatter(Ysel, Xsel, 20, 'r','filled');
        end
        drawnow;
        frame = getframe(gca);
        % Resize frame.cdata to [nRows nCols 3] to avoid size mismatches
        resizedFrame = imresize(frame.cdata, [nRows, nCols]);
        allOverlayImages(:,:,:,s) = resizedFrame;
        close(fig);
    end
    figMontage = figure('Visible','off');
    montage(allOverlayImages);
    montageFile = fullfile(overlayDir, 'GlobalAIF_AllSlices_Montage.png');
    exportgraphics(gca, montageFile, 'Resolution', 300);
    close(figMontage);
end
