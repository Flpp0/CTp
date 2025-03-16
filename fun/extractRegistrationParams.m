function [translations, rotations, translations_mm, time_values, highMovementsTable] = extractRegistrationParams(volumeTransforms, sortedPixelSpacing, uniqueTimePointsSeconds, upsample_factor)
% extractRegistrationParams Extracts translation and rotation parameters from volume transforms.
%
% Syntax:
%   [translations, rotations, translations_mm, time_values, highMovementsTable] = extractRegistrationParams(volumeTransforms, sortedPixelSpacing, uniqueTimePointsSeconds, upsample_factor)
%
% Description:
%   This function extracts translation and Euler angle rotation parameters from a cell array of 
%   rigidtform3d transformation objects. Translations are converted to millimeters using the pixel spacing.
%   It also identifies time points with high movement based on predefined
%   thresholds. High movement time points are elimanated from the further
%   steps. 
%
% Inputs:
%   volumeTransforms        - Cell array of rigidtform3d transformation objects.
%   sortedPixelSpacing      - Cell array of pixel spacing values (for the first slice).
%   uniqueTimePointsSeconds - Numeric array of unique time points (in seconds).
%   upsample_factor         - Upsampling factor (1 means no upsampling).
%
% Outputs:
%   translations    - minTimePoints x 3 array of translations.
%   rotations       - minTimePoints x 3 array of Euler angle rotations (in degrees).
%   translations_mm - Translations scaled to millimeters.
%   time_values     - Column vector of time indices.
%   highMovementsTable - Table summarizing indices and values for high movements.
%
% Example:
%   [trans, rot, trans_mm, time_vals, highMovements] = extractRegistrationParams(volumeTransforms, sortedPixelSpacing, uniqueTimePointsSeconds, 1);
%
% See also: rotationMatrixToEulerAngles

    minTimePoints = numel(volumeTransforms);
    translations = zeros(minTimePoints, 3);
    rotations = zeros(minTimePoints, 3);
    
    for t = 1:minTimePoints
        tform = volumeTransforms{t};
        R = tform.Rotation;
        T = tform.Translation;
        translations(t, :) = T;
        [rotX, rotY, rotZ] = rotationMatrixToEulerAngles(R);
        rotations(t, :) = [rotX, rotY, rotZ];
    end
    
    row_spacing = sortedPixelSpacing{1}(1);
    col_spacing = sortedPixelSpacing{1}(2);
    if numel(sortedPixelSpacing) > 1
        original_slice_spacing = abs(uniqueTimePointsSeconds(2) - uniqueTimePointsSeconds(1));
    else
        original_slice_spacing = 1;
    end
    new_slice_spacing = original_slice_spacing / upsample_factor;
    translations_mm = translations .* [col_spacing, row_spacing, new_slice_spacing];
    time_values = (1:minTimePoints).';
    
    translation_magnitudes = vecnorm(translations_mm, 2, 2);
    rotation_magnitudes = vecnorm(rotations, 2, 2);
    
    highTranslationThreshold = 5; % mm
    highRotationThreshold = 5;    % degrees
    highTranslationIndices = find(translation_magnitudes > highTranslationThreshold);
    highRotationIndices = find(rotation_magnitudes > highRotationThreshold);
    highMovementIndices = unique([highTranslationIndices; highRotationIndices]);
    
    highMovementsTable = table(highMovementIndices, time_values(highMovementIndices), ...
        translation_magnitudes(highMovementIndices), rotation_magnitudes(highMovementIndices), ...
        'VariableNames', {'Index', 'Time_s', 'Translation_mm', 'Rotation_deg'});
    
    disp('Registration parameters extracted and analyzed.');
end
