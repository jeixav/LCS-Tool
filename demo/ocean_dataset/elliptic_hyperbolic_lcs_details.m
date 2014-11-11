%% Input parameters
domain = [0,6;-34,-28];
timespan = [100,130];
resolutionX = 400;
resolutionY = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

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
lDerivative = @(t,x,~)derivative(t,x,vlon_interpolant,vlat_interpolant);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgEigenvalueFromMainGrid = false;
cgAuxGridRelDelta = .01;

% Lambda-lines
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6,'initialStep',1e-2);
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6];
poincareSection(2).endPosition = [1.3,-30.9;1.9,-31.1];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
% Set maximum orbit length to twice the expected circumference
for i = 1:nPoincareSection
    % FIXME Poincare section end point is not always close to vortex
    % centre, therefore rOrbit is not a good approximation of radius
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end

% Strainlines
strainlineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLocalMaxDistance = 2*gridSpace;
strainlineOdeSolverOptions = odeset('relTol',1e-6);

% Stretchlines
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 4*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-6);

% Graphics properties
strainlineColor = 'r';
stretchlineColor = 'b';
lambdaLineColor = [0,.6,0];
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
% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection,'UniformOutput',false);
hPoincareSection = [hPoincareSection{:}];
set(hPoincareSection,'color',lambdaLineColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineColor)
set(hPoincareSection,'MarkerEdgeColor','w')
drawnow

[etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',true);

% Plot all closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for i = 1:nPoincareSection
    hClosedLambdaLinePos{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{1},'UniformOutput',false);
    hClosedLambdaLineNeg{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{2},'UniformOutput',false);
end
hClosedLambdaLine = vertcat(vertcat(hClosedLambdaLinePos{:}),vertcat(hClosedLambdaLineNeg{:}));
hClosedLambdaLine = [hClosedLambdaLine{:}];
set(hClosedLambdaLine,'color',lambdaLineColor)

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2),'UniformOutput',false);
hLambdaLineLcsPos = [hLambdaLineLcsPos{:}];
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2),'UniformOutput',false);
hLambdaLineLcsNeg = [hLambdaLineLcsNeg{:}];
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineColor)
set(hLambdaLineLcs,'linewidth',2)
drawnow

%% Hyperbolic strainline LCSs
[strainlineLcs,strainlineLcsInitialPosition] = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineOdeSolverOptions);

% Remove strainlines inside elliptic regions
for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlineLcs = remove_strain_in_elliptic(strainlineLcs,closedLambdaLine{i}{1}{end});
    strainlineLcs = remove_strain_in_elliptic(strainlineLcs,closedLambdaLine{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(strainlineLcsInitialPosition(1,:),strainlineLcsInitialPosition(2,:),closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2));
    strainlineLcsInitialPosition = strainlineLcsInitialPosition(:,~idx);
    idx = inpolygon(strainlineLcsInitialPosition(1,:),strainlineLcsInitialPosition(2,:),closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2));
    strainlineLcsInitialPosition = strainlineLcsInitialPosition(:,~idx);
end

% Plot hyperbolic strainline LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs,'UniformOutput',false);
hStrainlineLcs = [hStrainlineLcs{:}];
set(hStrainlineLcs,'color',strainlineColor)
hStrainlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineLcsInitialPosition(1,idx),strainlineLcsInitialPosition(2,idx)),1:size(strainlineLcsInitialPosition,2),'UniformOutput',false);
hStrainlineLcsInitialPosition = [hStrainlineLcsInitialPosition{:}];
set(hStrainlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStrainlineLcsInitialPosition,'marker','o')
set(hStrainlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineLcsInitialPosition,'MarkerFaceColor',strainlineColor)

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
hPoincareSection = copyobj(hPoincareSection,hAxes);

% Plot all closed lambda lines
hClosedLambdaLine = copyobj(hClosedLambdaLine,hAxes);

% Plot lambda-line LCSs
hLambdaLineLcs = copyobj(hLambdaLineLcs,hAxes);
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions);

% Remove stretchlines inside elliptic regions
for i = 1:nPoincareSection
    % Remove stretchlines inside elliptic regions
    stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,closedLambdaLine{i}{1}{end});
    stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,closedLambdaLine{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(stretchlineLcsInitialPosition(1,:),stretchlineLcsInitialPosition(2,:),closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2));
    stretchlineLcsInitialPosition = stretchlineLcsInitialPosition(:,~idx);
    idx = inpolygon(stretchlineLcsInitialPosition(1,:),stretchlineLcsInitialPosition(2,:),closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2));
    stretchlineLcsInitialPosition = stretchlineLcsInitialPosition(:,~idx);
end

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs,'UniformOutput',false);
hStretchlineLcs = [hStretchlineLcs{:}];
set(hStretchlineLcs,'color',stretchlineColor)
hStretchlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineLcsInitialPosition(1,idx),stretchlineLcsInitialPosition(2,idx)),1:size(stretchlineLcsInitialPosition,2),'UniformOutput',false);
hStretchlineLcsInitialPosition = [hStretchlineLcsInitialPosition{:}];
set(hStretchlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStretchlineLcsInitialPosition,'marker','o')
set(hStretchlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchlineLcsInitialPosition,'MarkerFaceColor',stretchlineColor)

uistack(hLambdaLineLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
