% Some LCS Tool calculations require that the resolution in X and Y be
% equal. This function calculates the y-resolution that gives grid point
% spacing that is as close as possible to grid point spacing in the
% x-direction

function [resolutionY,deltaX] = equal_resolution(domain,resolutionX)

deltaX = (domain(1,2) - domain(1,1))/(resolutionX - 1);

resolutionY = round((domain(2,2) - domain(2,1))/deltaX) + 1;
