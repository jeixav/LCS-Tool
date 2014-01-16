%% Input parameters
domain = [0,6;-34,-28];
resolution = [400,400];
timespan = [98,128];

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
cgEigenvalueFromMainGrid = true;
cgAuxGridRelDelta = 0.1;

% Lambda-lines
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);

% Strainlines
strainlineLcsMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLcsLocalMaxDistance = 2*gridSpace;
strainlineLcsOdeSolverOptions = odeset('relTol',1e-4);

% Stretchlines
stretchlineLcsMaxLength = 20;
stretchlineLcsLocalMaxDistance = 4*gridSpace;
stretchlineLcsOdeSolverOptions = odeset('relTol',1e-4);

% Graphics properties
strainlineLcsColor = 'r';
stretchlineLcsColor = 'b';
lambdaLineLcsColor = [0,.6,0];
lcsInitialPositionMarkerSize = 2;

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
ftle_(isnan(ftle_)) = max(ftle_(:));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
drawnow

%% Lambda-line LCSs
% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [1.5,-30.9;1.9,-31.1];
poincareSection(2).endPosition = [2.9,-29.2;3.2,-29];

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',lambdaLineLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')
drawnow

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    % FIXME Poincare section end point is not always close to vortex
    % centre, therefore rOrbit is not a good approximation of radius
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end

[shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
s = warning('off','integrate_line:isDiscontinuousLargeAngle');
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',true);
warning(s)

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineLcsColor)
set(hLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaLinePos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{1});
    hClosedLambdaLineNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{2});
end
hClosedLambdaLine = horzcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
set(hClosedLambdaLine,'color',lambdaLineLcsColor)
drawnow

%% Hyperbolic strainline LCSs
[strainlineLcs,strainlineLcsInitialPosition] = seed_curves_from_lambda_max(strainlineLcsLocalMaxDistance,strainlineLcsMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineLcsOdeSolverOptions);

% Plot hyperbolic strainline LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs);
set(hStrainlineLcs,'color',strainlineLcsColor)
hStrainlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineLcsInitialPosition(1,idx),strainlineLcsInitialPosition(2,idx)),1:size(strainlineLcsInitialPosition,2));
set(hStrainlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStrainlineLcsInitialPosition,'marker','o')
set(hStrainlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineLcsInitialPosition,'MarkerFaceColor',strainlineLcsColor)

uistack(hLambdaLineLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(domain);
title(hAxes,'Stretchline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

% Plot finite-time Lyapunov exponent
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',lambdaLineLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineLcsColor)
set(hLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaLinePos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{1});
    hClosedLambdaLineNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{2});
end
hClosedLambdaLine = horzcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
set(hClosedLambdaLine,'color',lambdaLineLcsColor)
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(stretchlineLcsLocalMaxDistance,stretchlineLcsMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineLcsOdeSolverOptions);

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs);
set(hStretchlineLcs,'color',stretchlineLcsColor)
hStretchlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineLcsInitialPosition(1,idx),stretchlineLcsInitialPosition(2,idx)),1:size(stretchlineLcsInitialPosition,2));
set(hStretchlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStretchlineLcsInitialPosition,'marker','o')
set(hStretchlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchlineLcsInitialPosition,'MarkerFaceColor',stretchlineLcsColor)

uistack(hLambdaLineLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
