function smoothedImages = applyGaussianFiltering(imagesCellArray, sigma)
% applyGaussianFiltering Applies a Gaussian filter to each image in the input cell array.
%
% Syntax:
%   smoothedImages = applyGaussianFiltering(imagesCellArray, sigma)
%
% Description:
%   This function iterates through a cell array of images and applies a Gaussian filter
%   with a specified sigma to each non-empty image.
%
% Inputs:
%   imagesCellArray - Cell array of images.
%   sigma           - Standard deviation for the Gaussian kernel.
%
% Output:
%   smoothedImages - Cell array of filtered images.
%
% Example:
%   filteredImages = applyGaussianFiltering(myImages, 2);
%
% See also: imgaussfilt

    numSlices = size(imagesCellArray, 1);
    numTimePoints = size(imagesCellArray, 2);
    smoothedImages = cell(size(imagesCellArray));

    for s = 1:numSlices
        for t = 1:numTimePoints
            if ~isempty(imagesCellArray{s, t})
                smoothedImages{s, t} = imgaussfilt(imagesCellArray{s, t}, sigma);
            else
                smoothedImages{s, t} = [];
            end
        end
    end
end
