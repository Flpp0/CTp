function vpctDescription = getVPCTDescription(uniqueSeriesDescriptions)
% getVPCTDescription Retrieves the VPCT series description from a list.
%
% Syntax:
%   vpctDescription = getVPCTDescription(uniqueSeriesDescriptions)
%
% Description:
%   Searches through a cell array of series description strings and returns
%   the first string that contains the substring 'VPCT' (case-insensitive).
%
% Inputs:
%   uniqueSeriesDescriptions - Cell array of series description strings.
%
% Output:
%   vpctDescription - String containing the VPCT series description.
%
% Example:
%   vpctDescription = getVPCTDescription({'CT', 'VPCT Perfusion 5.0', 'MR'});


    vpctDescription = '';
    for i = 1:length(uniqueSeriesDescriptions)
        if contains(uniqueSeriesDescriptions{i}, 'VPCT', 'IgnoreCase', true)
            vpctDescription = uniqueSeriesDescriptions{i};
            break;
        end
    end
    if isempty(vpctDescription)
        error('No VPCT series description found.');
    end
end
