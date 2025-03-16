function downsampledVolume = downsampleXY(volume, factor)
% downsampleXY Downsamples a 3D volume in the XY plane only.
%
% Syntax:
%   downsampledVolume = downsampleXY(volume, factor)
%
% Description:
%   This function downsamples a 3D volume by the specified factor in the XY plane while 
%   leaving the Z-dimension (slices) unchanged. Bilinear interpolation is used for resizing.
%
% Inputs:
%   volume - A 3D array representing the input volume.
%   factor - Downsampling factor (e.g., a factor of 2 reduces each XY dimension by half).
%
% Output:
%   downsampledVolume - The downsampled 3D volume.

    [rows, cols, slices] = size(volume);
    newRows = floor(rows / factor);
    newCols = floor(cols / factor);
    downsampledVolume = zeros(newRows, newCols, slices, 'like', volume);
    for z = 1:slices
        downsampledVolume(:, :, z) = imresize(volume(:, :, z), 1 / factor, 'bilinear');
    end
end
