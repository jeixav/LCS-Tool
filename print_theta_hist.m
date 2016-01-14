% print_theta_hist Print a list of angular change values
%
% SYNTAX
% print_theta_hist(theta)
%
% EXAMPLE
% epsilon = .1;
% amplitude = .1;
% omega = pi/5;
% domain = [0,2;0,1];
% resolution = [750,375];
% timespan = [0,5];
% addpath(fullfile('demo','double_gyre'))
% lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
% incompressible = true;
% [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);
% xi1 = cgEigenvector(:,[1,2]);
% xi1 = shiftdim(reshape(xi1,[fliplr(resolution),2]),2);
% [thetaX,thetaY] = angle_change(xi1);
% print_theta_hist([thetaX(:);thetaY(:)])

function print_theta_hist(theta)

[n,xout] = hist(theta,deg2rad(5:10:175));

for i = 1:length(n)
    fprintf('%g\t%g\t%g\n',rad2deg(xout(i)),n(i),n(i)/sum(n))
end
