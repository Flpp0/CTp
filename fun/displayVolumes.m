function displayVolumes(smoothedHUImages, registeredHUImages, uniqueSliceLocations)
% displayVolumes Displays original and registered volumes for selected slices.
%
% Syntax:
%   displayVolumes(smoothedHUImages, registeredHUImages, uniqueSliceLocations)
%
% Description:
%   This function displays side-by-side the original (smoothed) and registered volumes for each slice
%   using sliceViewer. It constructs 3D volumes from the cell arrays and then displays them.
%
% Inputs:
%   smoothedHUImages - Cell array of original images (nSlices x nTimePoints).
%   registeredHUImages - Cell array of registered images (nSlices x nTimePoints).
%   uniqueSliceLocations - Array of unique slice location values.
%
% Example:
%   displayVolumes(smoothedHUImages, registeredHUImages, uniqueSliceLocations);


    numSlices = numel(uniqueSliceLocations);
    registered_volumes = cell(numSlices,1);
    original_volumes = cell(numSlices,1);
    slices_to_compare = 1:numSlices;
    compareSlices = true;
    
    if compareSlices
        for slice_idx = 1:numSlices
            current_slice = uniqueSliceLocations(slice_idx);
            valid_time_points = find(~cellfun('isempty', registeredHUImages(slice_idx, :)));
            num_valid_times = numel(valid_time_points);
            if num_valid_times == 0
                warning('No valid registered images for slice %.2f. Skipping visualization.', current_slice);
                continue;
            end
            [rows, cols] = size(registeredHUImages{slice_idx, valid_time_points(1)});
            slice_location_vol_reg = zeros(rows, cols, num_valid_times);
            slice_location_vol_ori = zeros(rows, cols, num_valid_times);
            for t = 1:num_valid_times
                time_idx = valid_time_points(t);
                slice_location_vol_reg(:,:,t) = registeredHUImages{slice_idx, time_idx};
                slice_location_vol_ori(:,:,t) = smoothedHUImages{slice_idx, time_idx};
            end
            registered_volumes{slice_idx} = slice_location_vol_reg;
            original_volumes{slice_idx} = slice_location_vol_ori;
            if ismember(slice_idx, slices_to_compare)
                fig1 = figure('Name', sprintf('Original Volume - Slice %.2f', current_slice));
                sliceViewer(slice_location_vol_ori);
                title(sprintf('Original Volume - Slice %.2f', current_slice));
                colormap('gray'); colorbar;
                fig2 = figure('Name', sprintf('Registered Volume - Slice %.2f', current_slice));
                sliceViewer(slice_location_vol_reg);
                title(sprintf('Registered Volume - Slice %.2f', current_slice));
                colormap('gray'); colorbar;
            end
        end
    end
end
