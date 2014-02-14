%% Input parameters
domain = [0,6;-34,-28];
resolution = [400,400];
timespan = [100,130];

%% Velocity definition
load('ocean_geostrophic_velocity.mat')
% Set velocity to zero at boundaries
vlon(:,[1,end],:) = 0;
vlon(:,:,[1,end]) = 0;
vlat(:,[1,end],:) = 0;
vlat(:,:,[1,end]) = 0;
interpMethod = 'spline';
vlon_interpolant = griddedInterpolant({time,lat,lon},vlon,interpMethod);
vlat_interpolant = griddedInterpolant({time,lat,lon},vlat,interpMethod);
lDerivative = @(t,x,~)flowdata_derivative(t,x,vlon_interpolant,vlat_interpolant);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgEigenvalueFromMainGrid = false;
cgAuxGridRelDelta = 0.01;

% Lambda-lines
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6,'initialStep',1e-2);

% Strainlines
strainlineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLocalMaxDistance = 2*gridSpace;
strainlineOdeSolverOptions = odeset('relTol',1e-4);

% Stretchlines
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 4*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-4);

% Graphics properties
strainlineColor = 'r';
stretchlineColor = 'b';
lambdaLineColor = [0,.6,0];

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

%% Lambda-line LCSs
% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6];
poincareSection(2).endPosition = [1.3,-30.9;1.9,-31.1];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end

[shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions);

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineColor)
set(hLambdaLineLcs,'linewidth',2)
drawnow

%% Hyperbolic strainline LCSs
strainlineLcs = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineOdeSolverOptions);

% Remove strainlines inside elliptic regions
for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedLambdaLine{i}{1}{1});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedLambdaLine{i}{2}{1});   
    % Remove initial positions inside elliptic regions
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
end

% Plot hyperbolic strainline LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs);
set(hStrainlineLcs,'color',strainlineColor)
uistack(hLambdaLineLcs,'top')
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(domain);
title(hAxes,'Stretchline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineColor)
set(hLambdaLineLcs,'linewidth',2)
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
stretchlineLcs = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions);

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs);
set(hStretchlineLcs,'color',stretchlineColor)

uistack(hLambdaLineLcs,'top')
