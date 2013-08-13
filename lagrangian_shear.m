% lagrangian_shear Lagrangian shear vector
%
% SYNTAX
% [etaPos,etaNeg] = lagrangian_shear(cgEigenvector,cgEigenvalue)
%
% DESCRIPTION
% Calculate the Lagragian shear vectors as defined in Equation 18 of
% DOI:10.1016/j.physd.2012.06.012

function [etaPos,etaNeg] = lagrangian_shear(cgEigenvector,cgEigenvalue)

l1 = cgEigenvalue(:,1);
l2 = cgEigenvalue(:,2);
xi1 = cgEigenvector(:,1:2);
xi2 = cgEigenvector(:,3:4);

etaPos(:,1) = sqrt(sqrt(l2)./(sqrt(l1) + sqrt(l2))).*xi1(:,1) + sqrt(sqrt(l1)./(sqrt(l1) + sqrt(l2))).*xi2(:,1);
etaPos(:,2) = sqrt(sqrt(l2)./(sqrt(l1) + sqrt(l2))).*xi1(:,2) + sqrt(sqrt(l1)./(sqrt(l1) + sqrt(l2))).*xi2(:,2);

etaNeg(:,1) = sqrt(sqrt(l2)./(sqrt(l1) + sqrt(l2))).*xi1(:,1) - sqrt(sqrt(l1)./(sqrt(l1) + sqrt(l2))).*xi2(:,1);
etaNeg(:,2) = sqrt(sqrt(l2)./(sqrt(l1) + sqrt(l2))).*xi1(:,2) - sqrt(sqrt(l1)./(sqrt(l1) + sqrt(l2))).*xi2(:,2);
