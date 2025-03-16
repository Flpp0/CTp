function vofStruct = computeGlobalVOF(allSlicesData, pAUC, pIreg, k, outDir, common_time_base, meanSlices, brainMasks, axes, centerOfMass)
% computeGlobalVOF Computes the global Venous Output Function (VOF) using a two-step filtering and clustering approach.
%
% Syntax:
%   vofStruct = computeGlobalVOF(allSlicesData, pAUC, pIreg, k, outDir, common_time_base, meanSlices, brainMasks, axes, centerOfMass)
%
% Description:
%   This function computes a single global VOF from the tissue concentration curves (TCCs)
%   stored in allSlicesData, but only for voxels located in the back of the brain (Quadrants 1 and 2).
%   The process includes:
%     1. Aggregating filtered TCCs and corresponding spatial coordinates from the slices.
%     2. Applying "second AUC" filtering (keeping only the top (100-pAUC)% of curves).
%     3. Performing roughness filtering by removing the top pIreg% of rough curves.
%     4. Conducting a two-step k-means clustering on the raw TACs curves to identify the cluster 
%        with the highest and latest peak.
%     5. Mapping the final cluster back to the original TACs amplitudes to compute the global VOF.
%     6. Generating and saving overlay visualizations per slice.
%
%   The VOF computation is based on a method similar to the approach for automatic selection of 
%   the arterial input function described in:
%
%       K. Mouridsen, S. Christensen, L. Gyldensted, and L. Østergaard, “Automatic selection of 
%       arterial input function using cluster analysis,” Magnetic Resonance in Med, vol. 55, no. 3, 
%       pp. 524–531, Mar. 2006, doi: 10.1002/mrm.20759.
%
% Inputs:
%   allSlicesData    - Cell array (one element per slice) with fields:
%                         .TCCs_filtered (V x T): Filtered time-attenuation curves.
%                         .voxIdx        (V x 1): Voxel indices within the slice.
%                         .rows, .cols   : Dimensions of the slice.
%                         .mask          : Binary mask (nRows x nCols).
%   pAUC             - Percentage threshold (0 to 100) for "second AUC" filtering.
%   pIreg            - Percentage threshold (0 to 100) for roughness filtering.
%   k                - Number of clusters for k-means clustering.
%   outDir           - Output directory to save diagnostic plots.
%   common_time_base - 1 x T vector of common time points.
%   meanSlices       - 3D array (nRows x nCols x numSlices) of mean precontrast images.
%   brainMasks       - 3D logical array (nRows x nCols x numSlices) representing brain masks.
%   axes             - Structure with principal axes (e.g., LeftRight, FrontBack).
%   centerOfMass     - [X, Y, Z] vector of the brain's center of mass.
%
% Output:
%   vofStruct - Structure containing:
%                .VOFTime   - The common time base.
%                .VOF       - The global VOF curve (raw TAC mean from the final cluster).
%                .voxelList - Global indices (in the concatenated TCC array) of the final cluster.
%                .Xcoord    - Row coordinates of the final cluster voxels.
%                .Ycoord    - Column coordinates of the final cluster voxels.
%                .sliceIdx  - Slice indices corresponding to each voxel.
%
% See also: computeRoughness, computeFirstMoment

    %% 1) Gather TCCs from All Slices (Only Back Voxels: Quadrants 1 and 2)
    C_original_all = [];
    Xcoord_all     = [];
    Ycoord_all     = [];
    sliceIdx_all   = [];
    bigIdx_all     = [];
    offset         = 0;
    
    numSlices = numel(allSlicesData);
    for sIdx = 2:numSlices-1  % Exclude first and last slices due to artifacts
        sliceData = allSlicesData{sIdx};
        if isempty(sliceData) || ~isfield(sliceData, 'TCCs_filtered') || isempty(sliceData.TCCs_filtered)
            continue;
        end
        
        % ----- Quadrant Filtering: Select voxels in Quadrants 1 and 2 -----
        backVoxelsMask = (sliceData.quadrant == 1) | (sliceData.quadrant == 2);
        if ~any(backVoxelsMask)
            continue;
        end
        
        TACs   = sliceData.TCCs_filtered(backVoxelsMask, :);
        voxIdx = sliceData.voxIdx(backVoxelsMask);
        rows   = sliceData.rows;
        cols   = sliceData.cols;
        quadrants = sliceData.quadrant(backVoxelsMask);
        
        if isempty(TACs)
            continue;
        end
        
        [Xc, Yc] = ind2sub([rows, cols], voxIdx);
        nV = size(TACs,1);
        C_original_all = [C_original_all; TACs];
        newGlobalIdx = (offset + 1 : offset + nV).';
        bigIdx_all = [bigIdx_all; newGlobalIdx];
        Xcoord_all = [Xcoord_all; Xc];
        Ycoord_all = [Ycoord_all; Yc];
        sliceIdx_tmp = repmat(sIdx, [nV, 1]);
        sliceIdx_all = [sliceIdx_all; sliceIdx_tmp];
        offset = offset + nV;
    end
    
    if isempty(C_original_all)
        warning('No TCCs found across all slices after quadrant filtering. Returning empty struct.');
        vofStruct = [];
        return;
    end
    
    C_original = C_original_all;
    Xcoord = Xcoord_all;
    Ycoord = Ycoord_all;
    sliceIdxV = sliceIdx_all;
    bigIdx = bigIdx_all;
    
    if size(C_original,2) ~= length(common_time_base)
        warning('C_original time dimension mismatch with common_time_base.');
    end
    
    %% 2) AUC Filtering
    AUC = trapz(common_time_base, C_original, 2);
    Nvox = size(C_original,1);
    [~, idxAUCsort] = sort(AUC, 'ascend');
    cutAUC = round(Nvox * (pAUC / 100));
    if cutAUC < Nvox
        idxKeepAUC = idxAUCsort(cutAUC+1:end);
    else
        warning('No curves remain after VOF AUC filtering...');
        vofStruct = [];
        return;
    end
    
    C = C_original(idxKeepAUC, :);
    Xcoord = Xcoord(idxKeepAUC);
    Ycoord = Ycoord(idxKeepAUC);
    sliceIdxV = sliceIdxV(idxKeepAUC);
    bigIdx = bigIdx(idxKeepAUC);
    AUC_sel = AUC(idxKeepAUC);
    
    if isempty(C)
        warning('No curves remain after AUC filtering...');
        vofStruct = [];
        return;
    end
    
    %% 3) Roughness Filtering
    areaScale = trapz(common_time_base, C, 2);
    Cunit_area = C ./ areaScale;
    
    roughVals = computeRoughness(Cunit_area, common_time_base);
    [~, idxRsort] = sort(roughVals, 'descend');
    cutR = round(length(idxRsort) * (pIreg / 100));
    
    if cutR < length(idxRsort)
        idxKeepR = idxRsort(cutR+1:end);
    else
        warning('No curves remain after roughness filtering for VOF...');
        vofStruct = [];
        return;
    end
    
    Cunit_area = Cunit_area(idxKeepR, :);
    C = C(idxKeepR, :);
    Xcoord = Xcoord(idxKeepR);
    Ycoord = Ycoord(idxKeepR);
    sliceIdxV = sliceIdxV(idxKeepR);
    bigIdx = bigIdx(idxKeepR);
    
    if isempty(C)
        warning('No curves remain after roughness filtering...');
        vofStruct = [];
        return;
    end
    
    %% 4) Two-Step k-means Clustering: Select Cluster with Highest and Latest Peak
    % First Clustering
    [clusterIdx1, Cmeans1] = kmeans(C, k, 'Replicates', 5, 'MaxIter', 1500);
    [maxPeak1, timeOfMax1] = getClusterPeakInfo(Cmeans1, common_time_base);
    maxIdx1 = selectCluster(maxPeak1, timeOfMax1);
    
    % Plot First Clustering Means
    cmap_full = autumn(256);
    figClust1 = figure('Visible','off');
    hold on;
    for kk = 1:k
        plot(common_time_base, Cmeans1(kk,:), 'LineWidth', 2);
    end
    plot(common_time_base, Cmeans1(maxIdx1,:), 'k', 'LineWidth', 3);
    title('First Clustering (Global VOF) - Mean Curves');
    xlabel('Time (s)'); ylabel('Raw Attenuation');
    hold off;
    saveas(figClust1, fullfile(outDir, 'GlobalVOF_FirstClustering_MeanCurves.png'));
    close(figClust1);
    
    idxSubset1 = find(clusterIdx1 == maxIdx1);
    Csubset1 = Cunit_area(idxSubset1, :);
    CrawSubset1 = C(idxSubset1, :);
    bigIdx_sub1 = bigIdx(idxSubset1);
    Xsub1 = Xcoord(idxSubset1);
    Ysub1 = Ycoord(idxSubset1);
    Ssub1 = sliceIdxV(idxSubset1);
    
    % Second Clustering on the first subset
    [clusterIdx2, Cmeans2] = kmeans(CrawSubset1, k, 'Replicates', 5, 'MaxIter', 1000);
    [maxPeak2, timeOfMax2] = getClusterPeakInfo(Cmeans2, common_time_base);
    maxIdx2 = selectCluster(maxPeak2, timeOfMax2);
    
    % Plot Second Clustering Means
    figClust2 = figure('Visible','off');
    hold on;
    for kk = 1:k
        plot(common_time_base, Cmeans2(kk,:), 'LineWidth', 2);
    end
    plot(common_time_base, Cmeans2(maxIdx2,:), 'k', 'LineWidth', 3);
    title('Second Clustering (Global VOF) - Mean Curves');
    xlabel('Time (s)'); ylabel('Raw Attenuation');
    hold off;
    saveas(figClust2, fullfile(outDir, 'GlobalVOF_SecondClustering_MeanCurves.png'));
    close(figClust2);
    
    idxFinal2 = find(clusterIdx2 == maxIdx2);
    finalIdx_all = bigIdx_sub1(idxFinal2);
    Xfinal = Xsub1(idxFinal2);
    Yfinal = Ysub1(idxFinal2);
    Sfinal = Ssub1(idxFinal2);
    
    % Retrieve raw DSC curves for final voxels from the original big array
    finalClusterOriginal = C_original_all(finalIdx_all, :);
    VOFraw_unscaled = mean(finalClusterOriginal, 1);
    VOFtime = common_time_base;
    
    %% 5) Build Output Structure
    vofStruct = struct();
    vofStruct.VOFTime = VOFtime;
    vofStruct.VOF = VOFraw_unscaled;
    vofStruct.voxelList = finalIdx_all;
    vofStruct.Xcoord = Xfinal;
    vofStruct.Ycoord = Yfinal;
    vofStruct.sliceIdx = Sfinal;
    
    fprintf('computeGlobalVOF: Final cluster has %d voxels.\n', numel(finalIdx_all));
    
    %% 6) Overlays Per Slice (Visualizing VOF Voxels)
    overlayDir = fullfile(outDir, 'GlobalVOF_VoxelOverlays');
    if ~exist(overlayDir, 'dir')
        mkdir(overlayDir);
    end
    
    [nRows, nCols, nSlcs] = size(meanSlices);
    for s = 1:nSlcs
        thisSliceMask = (Sfinal == s);
        if ~any(thisSliceMask)
            continue;
        end
        Xsel = Xfinal(thisSliceMask);
        Ysel = Yfinal(thisSliceMask);
        figName = sprintf('Slice_%d_GlobalVOF_voxels', s);
        fig = figure('Visible','off');
        imshow(meanSlices(:,:,s), []);
        hold on;
        % Changed color from 'r' to 'b' below:
        scatter(Ysel, Xsel, 20, 'b', 'filled');
        title(sprintf('Global VOF Voxels - Slice %d', s));
        hold off;
        outFile = fullfile(overlayDir, [figName, '.png']);
        saveas(fig, outFile);
        close(fig);
        fprintf('Saved overlay: %s\n', outFile);
    end

    %% 7) Save a montage of all slices with marked voxels for VOF
    allOverlayImages = zeros(nRows, nCols, 3, nSlcs, 'uint8');
    for s = 1:nSlcs
        fig = figure('Visible','off');
        imshow(meanSlices(:,:,s), []);
        hold on;
        thisSliceMask = (Sfinal == s);
        if any(thisSliceMask)
            Xsel = Xfinal(thisSliceMask);
            Ysel = Yfinal(thisSliceMask);
            scatter(Ysel, Xsel, 20, 'b','filled');
        end
        drawnow;
        frame = getframe(gca);
        % Resize frame to ensure it matches [nRows, nCols, 3]
        resizedFrame = imresize(frame.cdata, [nRows, nCols]);
        allOverlayImages(:,:,:,s) = resizedFrame;
        close(fig);
    end
    figMontage = figure('Visible','off');
    montage(allOverlayImages);
    montageFile = fullfile(overlayDir, 'GlobalVOF_AllSlices_Montage.png');
    exportgraphics(gca, montageFile, 'Resolution', 300);
    close(figMontage);

    %% Helper Function: getClusterPeakInfo
    function [maxPeak, timeOfMax] = getClusterPeakInfo(Cmeans, timeBase)
        % getClusterPeakInfo Computes the maximum peak and corresponding time for each cluster.
        %
        % Syntax:
        %   [maxPeak, timeOfMax] = getClusterPeakInfo(Cmeans, timeBase)
        %
        % Description:
        %   For each cluster mean curve in Cmeans, this function finds the maximum value (peak)
        %   and the corresponding time (using timeBase). Both outputs are vectors of length equal
        %   to the number of clusters.
        %
        % Inputs:
        %   Cmeans  - Matrix (k x T) where each row is a cluster mean curve.
        %   timeBase- 1 x T vector of time points.
        %
        % Outputs:
        %   maxPeak   - k x 1 vector of maximum peak values for each cluster.
        %   timeOfMax - k x 1 vector of times corresponding to the maximum peak of each cluster.
        
        [maxPeak, peakIdx] = max(Cmeans, [], 2);
        timeOfMax = timeBase(peakIdx);
    end

    %% Helper Function: selectCluster
    function selectedIdx = selectCluster(maxPeak, timeOfMax)
        % selectCluster Selects a cluster based on the highest peak and latest peak time.
        %
        % Syntax:
        %   selectedIdx = selectCluster(maxPeak, timeOfMax)
        %
        % Description:
        %   Among clusters with the highest peak, this function selects the one with the latest time
        %   at which the peak occurs.
        %
        % Inputs:
        %   maxPeak   - k x 1 vector of peak values for each cluster.
        %   timeOfMax - k x 1 vector of times corresponding to each cluster's peak.
        %
        % Output:
        %   selectedIdx - Index of the selected cluster.
        
        maxPeakVal = max(maxPeak);
        highestPeakIdx = find(maxPeak == maxPeakVal);
        [~, latestPeakIdx] = max(timeOfMax(highestPeakIdx));
        selectedIdx = highestPeakIdx(latestPeakIdx);
    end

end
