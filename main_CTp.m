%% main.m - CTp
% Rigoni Filippo

clc; clear; close all;

% Add to path
% Ensure the CTP_matlab folder is added to path files



%% Phase 1: Data Acquisition & Preprocessing

%% 0. Configuration
baseDir = ''; % base directory in which we have the folders with the patients data (multiple patients at this location)
patientCode = ''; % patient code. Selecting one from the available in the 'baseDir'
folderPath = fullfile(baseDir, patientCode, 'DICOM'); % if needed, change 'DICOM'. Folder containing the DICOM files inside the specific patient 'patientCode'
upsample_factor = 1; % 1 means no upsampling. README for explanation.

if ~exist(folderPath, 'dir')
    error('The specified DICOM folder does not exist: %s', folderPath);
end

%% 1. Load DICOM Data
[info, images, pixelSpacing, uniqueSeriesDescriptions, dicomFiles] = loadDicomData(folderPath);

%% 2. Get VPCT Series Description
vpctDescription = getVPCTDescription(uniqueSeriesDescriptions);
disp('The VPCT description is:');
disp(vpctDescription);

%% 3. Load VPCT Images and Associated Pixel Spacing
[vpctImages, timePoints, sliceLocations, vpctPixelSpacing] = loadVPCTImages(info, dicomFiles, vpctDescription, pixelSpacing);

%% 4. Process PixelSpacing values for VPCT images
vpctPixelSpacing = fillNaNPixels(vpctPixelSpacing);

%% 5. Sort VPCT Images by Slice and Time
[uniqueTimePoints, uniqueTimePointsSeconds, uniqueSliceLocations, sortedHUImages, sortedPixelSpacing] = sortImages(vpctImages, timePoints, sliceLocations, vpctPixelSpacing);
clear vpctImages;

%% 6. (Optional) Apply Gaussian Filtering
applyGaussianFilteringFlag = false;
if applyGaussianFilteringFlag
    smoothedHUImages = applyGaussianFiltering(sortedHUImages, 0.5);
else
    smoothedHUImages = sortedHUImages;
end
clear sortedHUImages;




%% Phase 2: Volume Construction & Registration

%% 7. Construct Volumes and Perform 3D Rigid Registration
[volumes, nonEmptyTimeIndices, minTimePoints] = constructVolumes(smoothedHUImages);
interpolatedVolumes = upsampleVolumes(volumes, upsample_factor, baseDir, patientCode);
[registeredHUImages, volumeTransforms] = registerVolumes(interpolatedVolumes, uniqueSliceLocations, sortedPixelSpacing, baseDir, patientCode, upsample_factor);
clear volumes interpolatedVolumes;

%% 8. Extract and Plot Registration Parameters Over Time
[translations, rotations, translations_mm, time_values, highMovementsTable] = extractRegistrationParams(volumeTransforms, sortedPixelSpacing, uniqueTimePointsSeconds, upsample_factor);
plotRegistrationParameters(translations_mm, rotations, time_values, highMovementsTable, baseDir, patientCode, upsample_factor);

%% 9. Correct Positioning of Images in registeredHUImages
registeredHUImages = correctImagePositioning(smoothedHUImages, registeredHUImages, uniqueTimePointsSeconds, upsample_factor);

%% 10. Display Original vs. Registered Volumes for Selected Slices
displayVolumes(smoothedHUImages, registeredHUImages, uniqueSliceLocations);




%% Phase 3: Post-Registration Processing & Feature Extraction

%% 11. Apply Post-Registration Gaussian Filtering
gaussianSigma = 3;
smoothedHUImages = applyGaussianFiltering(registeredHUImages, gaussianSigma);

%% 12. Skull Stripping & Mean Slices Computation
[nRows, nCols] = deal(512,512);
[meanSlices, brainMasks] = performSkullStripping(smoothedHUImages, baseDir, patientCode, upsample_factor, nRows, nCols, 100);

%% 13. Compute or Load Principal Axes
[axes, centerOfMass, coeff, latent] = computePrincipalAxes(meanSlices, brainMasks, baseDir, patientCode, upsample_factor);

%% 14. Calculate and Save Tissue Concentration Curves (TCCs) with Outlier Removal
[common_time_base, allSlicesData, lastValidTimeSeconds, brainMasks] = calculateTCCs(registeredHUImages, smoothedHUImages, brainMasks, meanSlices, uniqueTimePointsSeconds, baseDir, patientCode, upsample_factor, axes, centerOfMass);




%% Phase 4: Perfusion Analysis & Mapping

%% 15. Automatic Global AIF/VOF Detection and Scaling
[aifStruct, vofStruct, globalAIFScaledVOF] = computeGlobalAIF_VOF(allSlicesData, common_time_base, meanSlices, brainMasks, axes, centerOfMass, baseDir, patientCode, upsample_factor);

%% 16. Initialize AIF and Compute Residue Functions using sSVD, cSVD, and oSVD
computeResidueFunctions(allSlicesData, common_time_base, globalAIFScaledVOF, upsample_factor, baseDir, patientCode);

%% 17. Compute Perfusion Maps for sSVD, cSVD, and oSVD Methods
computePerfusionMapsAll(baseDir, patientCode, brainMasks, meanSlices, upsample_factor);
close all; % Better to close all the figures before proceeding with the saving.

%% 18. Save Perfusion Maps for sSVD, cSVD, and oSVD Methods
savePerfusionMapImages(baseDir, patientCode, upsample_factor, meanSlices);

disp('Processing pipeline complete.');
