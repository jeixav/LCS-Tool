% lambda_line Null-geodesics of generalized Green-Lagrange Lorentzian
% metric
%
% SYNTAX
% [etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda)
%
% DESCRIPTION
% Calculate the quantity defined in Equation 3.2 of
% DOI:10.1017/jfm.2013.391.

function [etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda)

validateattributes(cgEigenvector,{'double'},{'size',[NaN,4]})
validateattributes(cgEigenvalue,{'double'},{'size',[NaN,2]})
validateattributes(lambda,{'double'},{'scalar'})

l1 = cgEigenvalue(:,1);
l2 = cgEigenvalue(:,2);
xi1 = cgEigenvector(:,1:2);
xi2 = cgEigenvector(:,3:4);

etaPos(:,1) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,1) + sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,1);
etaPos(:,2) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,2) + sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,2);

etaNeg(:,1) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,1) - sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,1);
etaNeg(:,2) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,2) - sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,2);
