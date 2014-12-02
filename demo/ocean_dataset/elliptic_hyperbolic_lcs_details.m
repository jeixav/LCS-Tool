%% Input parameters
timespan = [100,130];
domain = [0,6;-34,-28];
resolutionX = 400;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
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

% Lambda lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6];
poincareSection(2).endPosition = [1.3,-30.9;1.9,-31.1];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6,'initialStep',1e-2);
showPoincareGraph = true;

% Shrink lines
shrinkLineMaxLength = 20;
shrinkLineLocalMaxDistance = 2*deltaX;
shrinkLineOdeSolverOptions = odeset('relTol',1e-6);

% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 4*deltaX;
stretchLineOdeSolverOptions = odeset('relTol',1e-6);

% Graphic properties
repellingColor = 'r';
attractingColor = 'b';
ellipticColor = [0,.6,0];
initialPositionMarkerSize = 2;

hAxes = setup_figure(domain);
title(hAxes,'Repelling and elliptic LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
drawnow

%% Elliptic LCSs
% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection,'UniformOutput',false);
hPoincareSection = [hPoincareSection{:}];
set(hPoincareSection,'color',ellipticColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',ellipticColor)
set(hPoincareSection,'MarkerEdgeColor','w')
hPoincareSectionText = arrayfun(@(idx)text(poincareSection(idx).endPosition(2,1),poincareSection(idx).endPosition(2,2),['# ',num2str(idx)']),1:nPoincareSection,'UniformOutput',false);
hPoincareSectionText = [hPoincareSectionText{:}];
set(hPoincareSectionText,'parent',hAxes)
set(hPoincareSectionText,'color',ellipticColor)
drawnow

[etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',showPoincareGraph);

% Plot closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for i = 1:nPoincareSection
    hClosedLambdaLinePos{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{1},'UniformOutput',false);
    hClosedLambdaLineNeg{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{2},'UniformOutput',false);
end
hClosedLambdaLine = vertcat(vertcat(hClosedLambdaLinePos{:}),vertcat(hClosedLambdaLineNeg{:}));
hClosedLambdaLine = [hClosedLambdaLine{:}];
set(hClosedLambdaLine,'color',ellipticColor)

% Plot elliptic LCSs
hEllipticLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2),'UniformOutput',false);
hEllipticLcsPos = [hEllipticLcsPos{:}];
hEllipticLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2),'UniformOutput',false);
hEllipticLcsNeg = [hEllipticLcsNeg{:}];
hEllipticLcs = [hEllipticLcsPos,hEllipticLcsNeg];
set(hEllipticLcs,'color',ellipticColor)
set(hEllipticLcs,'linewidth',2)
drawnow

%% Repelling LCSs
[shrinkLine,shrinkLineInitialPosition] = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions);

% Remove shrink lines inside elliptic LCSs
for i = 1:nPoincareSection
    shrinkLine = remove_strain_in_elliptic(shrinkLine,closedLambdaLine{i}{1}{end});
    shrinkLine = remove_strain_in_elliptic(shrinkLine,closedLambdaLine{i}{2}{end});
    idx = inpolygon(shrinkLineInitialPosition(1,:),shrinkLineInitialPosition(2,:),closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2));
    shrinkLineInitialPosition = shrinkLineInitialPosition(:,~idx);
    idx = inpolygon(shrinkLineInitialPosition(1,:),shrinkLineInitialPosition(2,:),closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2));
    shrinkLineInitialPosition = shrinkLineInitialPosition(:,~idx);
end

% Plot repelling LCSs
hRepellingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),shrinkLine,'UniformOutput',false);
hRepellingLcs = [hRepellingLcs{:}];
set(hRepellingLcs,'color',repellingColor)
hShrinkLineInitialPosition = arrayfun(@(idx)plot(hAxes,shrinkLineInitialPosition(1,idx),shrinkLineInitialPosition(2,idx)),1:size(shrinkLineInitialPosition,2),'UniformOutput',false);
hShrinkLineInitialPosition = [hShrinkLineInitialPosition{:}];
set(hShrinkLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
set(hShrinkLineInitialPosition,'marker','o')
set(hShrinkLineInitialPosition,'MarkerEdgeColor','w')
set(hShrinkLineInitialPosition,'MarkerFaceColor',repellingColor)

uistack(hEllipticLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
uistack(hPoincareSectionText,'top')
drawnow

%% Attracting LCSs
hAxes = setup_figure(domain);
title(hAxes,'Attracting and elliptic LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

% Plot finite-time Lyapunov exponent
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))

% Copy objects from repelling LCS plot
hPoincareSection = copyobj(hPoincareSection,hAxes);
hPoincareSectionText = copyobj(hPoincareSectionText,hAxes);
hClosedLambdaLine = copyobj(hClosedLambdaLine,hAxes);
hEllipticLcs = copyobj(hEllipticLcs,hAxes);
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
[stretchLine,stretchLineInitialPosition] = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,closedLambdaLine{i}{1}{end});
    stretchLine = remove_strain_in_elliptic(stretchLine,closedLambdaLine{i}{2}{end});
    idx = inpolygon(stretchLineInitialPosition(1,:),stretchLineInitialPosition(2,:),closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2));
    stretchLineInitialPosition = stretchLineInitialPosition(:,~idx);
    idx = inpolygon(stretchLineInitialPosition(1,:),stretchLineInitialPosition(2,:),closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2));
    stretchLineInitialPosition = stretchLineInitialPosition(:,~idx);
end

% Plot attracting LCSs
hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchLine,'UniformOutput',false);
hAttractingLcs = [hAttractingLcs{:}];
set(hAttractingLcs,'color',attractingColor)
hStretchLineInitialPosition = arrayfun(@(idx)plot(hAxes,stretchLineInitialPosition(1,idx),stretchLineInitialPosition(2,idx)),1:size(stretchLineInitialPosition,2),'UniformOutput',false);
hStretchLineInitialPosition = [hStretchLineInitialPosition{:}];
set(hStretchLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
set(hStretchLineInitialPosition,'marker','o')
set(hStretchLineInitialPosition,'MarkerEdgeColor','w')
set(hStretchLineInitialPosition,'MarkerFaceColor',attractingColor)

uistack(hEllipticLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
uistack(hPoincareSectionText,'top')
