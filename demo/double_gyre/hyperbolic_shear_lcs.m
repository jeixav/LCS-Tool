%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;

incompressible = true;
derivative = @(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega);
domain = [-.1,2.1;-.05,1.05];
timespan = [0,20];
resolution = [551,276];

lambda = 1;

hyperbolicLcsMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLcsLocalMaxDistance = 2*gridSpace;

stretchlineLcsLocalMaxDistance = 4*gridSpace;

strainlineLcsColor = 'r';
stretchlineLcsColor = 'b';
lambdaLineLcsColor = [0,.6,0];

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvalue,cgEigenvector] = eig_cgStrain(derivative,domain,timespan,resolution,'incompressible',incompressible);

%% Lambda-line LCSs
% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [.5,.6;.35,.5];
poincareSection(2).endPosition = [1.5,.4;1.7,.5];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(80);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

[shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedOrbits = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection);

% Plot lambda-line LCSs
% η₊ outermost closed lambda-lines
hShearLcsPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
% η₋ outermost closed lambda-lines
hShearLcsNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
hShearLcs = [hShearLcsPos,hShearLcsNeg];
set(hShearLcs,'color',lambdaLineLcsColor)
set(hShearLcs,'linewidth',2)
drawnow

%% Hyperbolic strainline LCSs
strainlinePosition = seed_curves_from_lambda_max(strainlineLcsLocalMaxDistance,hyperbolicLcsMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution);

% Plot hyperbolic strainline LCSs
hStrainLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainLcs,'color',strainlineLcsColor)

uistack(hShearLcs,'top')
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(domain);
title(hAxes,'Stretchline and \lambda-line LCSs')

% Plot lambda-line LCSs
% η₊ outermost closed lambda-lines
hShearLcsPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
% η₋ outermost closed lambda-lines
hShearLcsNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
hShearLcs = [hShearLcsPos,hShearLcsNeg];
set(hShearLcs,'color',lambdaLineLcsColor)
set(hShearLcs,'linewidth',2)
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
stretchlinePosition = seed_curves_from_lambda_max(stretchlineLcsLocalMaxDistance,hyperbolicLcsMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution);

% Plot hyperbolic stretchline LCSs
hStretchLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlinePosition);
set(hStretchLcs,'color',stretchlineLcsColor)

uistack(hShearLcs,'top')
