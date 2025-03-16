function val = computeRoughness(Cmat, timeVec)
% computeRoughness Computes the roughness of each curve by integrating the squared second derivative.
%
% Syntax:
%   val = computeRoughness(Cmat, timeVec)
%
% Description:
%   For each row in Cmat (representing a curve), this function computes the second derivative with 
%   respect to time (given by timeVec), squares it, and integrates it over time to yield a roughness value.
%   Roughness = ∫ [second derivative]^2 dt
%
% Inputs:
%   Cmat    - Matrix of curves (each row is a curve, already normalized by area).
%   timeVec - Vector of time points corresponding to the columns of Cmat.
%
% Output:
%   val - Column vector of roughness values for each curve.

    dCdt = gradient(Cmat, mean(diff(timeVec)), 2); % 1st derivative wrt time
    d2Cdt2 = gradient(dCdt, mean(diff(timeVec)), 2); % 2nd derivative
    val = trapz(timeVec, (d2Cdt2.^2), 2); % integrate wrt time
end