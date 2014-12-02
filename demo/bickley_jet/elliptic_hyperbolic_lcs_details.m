%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
timespan = [0,2*lengthX/u];
domain = [2e6,.5*lengthX;[-1,.25]*2.25*lengthY];
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

%% Velocity definition
perturbationCase = 3;
phiTimespan = [0,25];
phiInitial = [0,0];
phiSol = ode45(@d_phi,phiTimespan,phiInitial);
timeResolution = 1e5;
phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
phi1Max = max(phi1);
lDerivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-4);

% Lambda lines
poincareSection.endPosition = [6.5,-1.4;4.5,-3.5]*1e6;
[poincareSection.numPoints] = deal(100);
rOrbit = hypot(diff(poincareSection.endPosition(:,1)),diff(poincareSection.endPosition(:,2)));
poincareSection.orbitMaxLength = 2*(2*pi*rOrbit);
dThresh = 1e-3;
lambda = .995;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
showPoincareGraph = true;

% Shrink lines
shrinkLineMaxLength = 1e8;
shrinkLineLocalMaxDistance = 4*deltaX;
shrinkLineOdeSolverOptions = odeset('relTol',1e-4);

% Stretch lines
stretchLineMaxLength = 1e8;
stretchLineLocalMaxDistance = 8*deltaX;
stretchLineOdeSolverOptions = odeset('relTol',1e-4);

% Graphic properties
repellingColor = 'r';
attractingColor = 'b';
ellipticColor = [0,.6,0];
initialPositionMarkerSize = 2;

hAxes = setup_figure(domain);
title(hAxes,'Repelling and elliptic LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

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
drawnow

[etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'dThresh',dThresh,'showGraph',showPoincareGraph);

% Plot closed lambda lines
hClosedLambdaLinePos = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{1}{1},'UniformOutput',false);
hClosedLambdaLineNeg = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{1}{2},'UniformOutput',false);
hClosedLambdaLine = vertcat(hClosedLambdaLinePos{1},hClosedLambdaLineNeg{1});
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
shrinkLine = remove_strain_in_elliptic(shrinkLine,closedLambdaLine{1}{1}{end});
shrinkLine = remove_strain_in_elliptic(shrinkLine,closedLambdaLine{1}{2}{end});
idx = inpolygon(shrinkLineInitialPosition(1,:),shrinkLineInitialPosition(2,:),closedLambdaLine{1}{1}{end}(:,1),closedLambdaLine{1}{1}{end}(:,2));
shrinkLineInitialPosition = shrinkLineInitialPosition(:,~idx);
idx = inpolygon(shrinkLineInitialPosition(1,:),shrinkLineInitialPosition(2,:),closedLambdaLine{1}{2}{end}(:,1),closedLambdaLine{1}{2}{end}(:,2));
shrinkLineInitialPosition = shrinkLineInitialPosition(:,~idx);

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
drawnow

%% Attracting LCSs
hAxes = setup_figure(domain);
title(hAxes,'Attracting and elliptic LCSs')

% Plot finite-time Lyapunov exponent
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))

% Copy objects from repelling LCS plot
hPoincareSection = copyobj(hPoincareSection,hAxes);
hClosedLambdaLine = copyobj(hClosedLambdaLine,hAxes);
hEllipticLcs = copyobj(hEllipticLcs,hAxes);
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
[stretchLine,stretchLineInitialPosition] = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
stretchLine = remove_strain_in_elliptic(stretchLine,closedLambdaLine{1}{1}{end});
stretchLine = remove_strain_in_elliptic(stretchLine,closedLambdaLine{1}{2}{end});
idx = inpolygon(stretchLineInitialPosition(1,:),stretchLineInitialPosition(2,:),closedLambdaLine{1}{1}{end}(:,1),closedLambdaLine{1}{1}{end}(:,2));
stretchLineInitialPosition = stretchLineInitialPosition(:,~idx);
idx = inpolygon(stretchLineInitialPosition(1,:),stretchLineInitialPosition(2,:),closedLambdaLine{1}{2}{end}(:,1),closedLambdaLine{1}{2}{end}(:,2));
stretchLineInitialPosition = stretchLineInitialPosition(:,~idx);

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
