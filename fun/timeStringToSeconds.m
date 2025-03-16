function seconds = timeStringToSeconds(timeStr)
% timeStringToSeconds Converts a DICOM acquisition time string to seconds.
%
% Syntax:
%   seconds = timeStringToSeconds(timeStr)
%
% Description:
%   This function converts a time string in the HHMMSS.FFFFFF format (common in DICOM
%   metadata) into a total number of seconds.
%
% Inputs:
%   timeStr - A string representing the time in HHMMSS.FFFFFF format.
%
% Outputs:
%   seconds - Total time in seconds as a double.
%
% Example:
%   s = timeStringToSeconds('123456.789');

    % If the time string contains a fractional part, separate it
    if contains(timeStr, '.')
        [timeStrMain, fracStr] = strtok(timeStr, '.');
        fracSeconds = str2double(['0.' fracStr(2:end)]);
    else
        timeStrMain = timeStr;
        fracSeconds = 0;
    end

    % Pad the main time string to ensure it has at least 6 characters (HHMMSS)
    timeStrMain = pad(timeStrMain, 6, 'left', '0');

    % Extract hours, minutes, and seconds
    hours = str2double(timeStrMain(1:2));
    minutes = str2double(timeStrMain(3:4));
    seconds_part = str2double(timeStrMain(5:6));

    % Compute total seconds
    seconds = hours * 3600 + minutes * 60 + seconds_part + fracSeconds;
end
