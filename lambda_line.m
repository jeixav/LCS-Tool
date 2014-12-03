% lambda_line Null-geodesics of generalized Green-Lagrange Lorentzian
% metric
%
% SYNTAX
% [etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda)
% [etaPos,etaNeg] = lambda_line(...,'forceComplexNaN',forceComplexNaN)
%
% INPUT ARGUMENTS
% forceComplexNaN: logical to control whether complex etaPos and etaNeg are
% replaced by NaN
%
% DESCRIPTION
% Calculate the quantity defined in Equation 3.2 of
% DOI:10.1017/jfm.2013.391.

function [etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda,varargin)

p = inputParser;

% FIXME Make validationFcn common with eig_cgStrain
addRequired(p,'cgEigenvector',@(i)validateattributes(i,{'double'},{'size',[NaN,4]}))
addRequired(p,'cgEigenvalue',@(i)validateattributes(i,{'double'},{'size',[NaN,2]}))
addRequired(p,'lambda',@(i)validateattributes(i,{'double'},{'scalar'}))
addParameter(p,'forceComplexNaN',false,@(i)validateattributes(i,{'logical'},{'scalar'}))

parse(p,cgEigenvector,cgEigenvalue,lambda,varargin{:})

cgEigenvector = p.Results.cgEigenvector;
cgEigenvalue = p.Results.cgEigenvalue;
lambda = p.Results.lambda;
forceComplexNaN = p.Results.forceComplexNaN;

l1 = cgEigenvalue(:,1);
l2 = cgEigenvalue(:,2);
xi1 = cgEigenvector(:,1:2);
xi2 = cgEigenvector(:,3:4);

etaPos(:,1) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,1) + sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,1);
etaPos(:,2) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,2) + sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,2);

etaNeg(:,1) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,1) - sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,1);
etaNeg(:,2) = sqrt((l2 - lambda^2)./(l2 - l1)).*xi1(:,2) - sqrt((lambda^2 - l1)./(l2 - l1)).*xi2(:,2);

[r_,~] = find(imag(etaPos) | imag(etaNeg));
if r_
    warning([mfilename,':complexEta'],'%d etaPos and etaNeg elements with imaginary parts. Maximum imaginary part: %g',numel(unique(r_)),max(imag([etaPos(:);etaNeg(:)])))
    if forceComplexNaN
        etaPos(r_,:) = nan;
        etaNeg(r_,:) = nan;
    end
end
