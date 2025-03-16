function [rx, ry, rz] = rotationMatrixToEulerAngles(R)
% rotationMatrixToEulerAngles Converts a 3x3 rotation matrix to Euler angles (in degrees).
%
% Syntax:
%   [rx, ry, rz] = rotationMatrixToEulerAngles(R)
%
% Description:
%   This function converts a 3x3 rotation matrix R into Euler angles (rx, ry, rz) measured in degrees.
%   The computation is based on the standard formulas. When the computed value 'sy' (used to detect a
%   singularity) is very small, a special case is applied to avoid numerical instability.
%
% Input:
%   R - A 3x3 rotation matrix.
%
% Outputs:
%   rx - Euler angle corresponding to rotation about the X-axis (in degrees).
%   ry - Euler angle corresponding to rotation about the Y-axis (in degrees).
%   rz - Euler angle corresponding to rotation about the Z-axis (in degrees).
%
% Example:
%   R = [0 -1 0; 1 0 0; 0 0 1];
%   [rx, ry, rz] = rotationMatrixToEulerAngles(R);
%
% See also: rad2deg, atan2

    sy = sqrt(R(1,1)^2 + R(2,1)^2);
    singular = sy < 1e-6;

    if ~singular
        rx = atan2(R(3,2), R(3,3));
        ry = atan2(-R(3,1), sy);
        rz = atan2(R(2,1), R(1,1));
    else
        rx = atan2(-R(2,3), R(2,2));
        ry = atan2(-R(3,1), sy);
        rz = 0;
    end

    % Convert angles from radians to degrees
    rx = rad2deg(rx);
    ry = rad2deg(ry);
    rz = rad2deg(rz);
end
