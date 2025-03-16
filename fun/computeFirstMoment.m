function val = computeFirstMoment(Cmat, timeVec)
% computeFirstMoment Computes the first moment (centroid) of each curve.
%
% Syntax:
%   val = computeFirstMoment(Cmat, timeVec)
%
% Description:
%   This function computes the first moment of each row in Cmat, defined as the ratio of the weighted
%   integral of time (weighted by the curve) to the integral of the curve. It is used to capture the 
%   "center of mass" of each curve.
%   first moment = (∫ t*C(t) dt) / (∫ C(t) dt)
%
% Inputs:
%   Cmat    - Matrix where each row represents a curve.
%   timeVec - Vector of time points corresponding to the columns of Cmat.
%
% Output:
%   val - Column vector of first moment values for each curve.

    tC = bsxfun(@times, Cmat, timeVec);
    num = trapz(timeVec, tC, 2); % ∫ t*C(t) dt
    den = trapz(timeVec, Cmat, 2); % ∫ C(t) dt
    val = num ./ den;
end