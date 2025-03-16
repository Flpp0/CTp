function normalizedVolume = normalizeVolumeGlobal(volume, globalMin, globalMax)
% normalizeVolumeGlobal Normalizes a volume using global intensity values.
%
% Syntax:
%   normalizedVolume = normalizeVolumeGlobal(volume, globalMin, globalMax)
%
% Description:
%   This function normalizes the input 3D volume so that its intensity values are scaled
%   to the range [0, 1]. The normalization is done using the provided global minimum and maximum
%   intensity values.
%
% Inputs:
%   volume    - 3D array representing the volume to be normalized.
%   globalMin - The global minimum intensity value across all volumes.
%   globalMax - The global maximum intensity value across all volumes.
%
% Output:
%   normalizedVolume - The normalized volume, with intensity values in the range [0, 1].

    normalizedVolume = (volume - globalMin) / (globalMax - globalMin);
end