%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];

flow.imposeIncompressibility = true;
perturbationCase = 3;
flow.derivative = @(t,x,useEoV)derivative(t,x,useEoV,u,lengthX,lengthY,epsilon,perturbationCase);
flow.resolution = [500,200];
% magicNumber gives the domain an aspect ratio similar to that used in
% doi:10.1016/j.physd.2012.06.012 and ensures grid spacing is equal in the
% x and y directions. In doi:10.1016/j.physd.2012.06.012, 
% magicNumber = 2.2599.
magicNumber = .5*pi*earthRadius/lengthY*double(flow.resolution(2)-1)/double(flow.resolution(1)-1);
flow.domain = [0,lengthX;[-1,1]*magicNumber*lengthY];
flow.timespan = [0,4*lengthX/u];
flow.periodicBc = [true,false];

lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-4);

hyperbolicLcsMaxLength = 1e8;
strainlineOdeSolverOptions = odeset('relTol',1e-4);
gridSpace = diff(flow.domain(1,:))/(double(flow.resolution(1))-1);
localMaxDistance = 8*gridSpace;

strainlineLcsColor = 'r';
stretchlineLcsColor = 'b';
lambdaLineLcsColor = [0,.6,0];

hAxes = setup_figure(flow.domain);
title(hAxes,'Strainline and \lambda-line LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvalue,cgEigenvector] = eig_cgStrain(flow);

%% Lambda-line LCSs
% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [6.5,-1.4;5.2,-1.4]*1e6;
poincareSection(2).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection(3).endPosition = [3.25,1.5;1.4,2.6]*1e6;
poincareSection(4).endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection(5).endPosition = [1.65e7,1.5e6;1.5e7,2.6e6];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(80);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

[shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedOrbits = poincare_closed_orbit_multi(flow,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions);

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
strainlinePosition = seed_curves_from_lambda_max(localMaxDistance,hyperbolicLcsMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),flow.domain,flow.resolution,'periodicBc',flow.periodicBc,'odeSolverOptions',strainlineOdeSolverOptions);

% Plot hyperbolic strainline LCSs
hStrainLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainLcs,'color',strainlineLcsColor)

uistack(hShearLcs,'top')
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(flow.domain);
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
stretchlinePosition = seed_curves_from_lambda_max(localMaxDistance,hyperbolicLcsMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),flow.domain,flow.resolution);

% Plot hyperbolic stretchline LCSs
hStretchLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlinePosition);
set(hStretchLcs,'color',stretchlineLcsColor)

uistack(hShearLcs,'top')
