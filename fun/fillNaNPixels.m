function vpctPixelSpacing = fillNaNPixels(vpctPixelSpacing)
% fillNaNPixels Replaces NaN values in PixelSpacing with the mode of available spacings.
%
% Syntax:
%   vpctPixelSpacing = fillNaNPixels(vpctPixelSpacing)
%
% Description:
%   This function scans the PixelSpacing values in the cell array and replaces any
%   NaN values with the mode of the available row and column spacing values.
%
% Input:
%   vpctPixelSpacing - Cell array where each cell contains a two-element vector [row_spacing, col_spacing].
%
% Output:
%   vpctPixelSpacing - Updated cell array with NaN values replaced.
%
% Example:
%   vpctPixelSpacing = fillNaNPixels(vpctPixelSpacing);


    allRowSpacings = cellfun(@(x) x(1), vpctPixelSpacing);
    allColSpacings = cellfun(@(x) x(2), vpctPixelSpacing);
    modeRowSpacing = mode(allRowSpacings(~isnan(allRowSpacings)));
    modeColSpacing = mode(allColSpacings(~isnan(allColSpacings)));
    
    for i = 1:length(vpctPixelSpacing)
        if any(isnan(vpctPixelSpacing{i}))
            vpctPixelSpacing{i}(isnan(vpctPixelSpacing{i})) = [modeRowSpacing, modeColSpacing];
            fprintf('Replaced NaN PixelSpacing for VPCT file %d with [%f, %f].\n', i, modeRowSpacing, modeColSpacing);
        end
    end
end
