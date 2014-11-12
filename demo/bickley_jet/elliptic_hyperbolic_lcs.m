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

% Lambda-lines
lambda = .995;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
poincareSection.endPosition = [6.5,-1.4;4.5,-3.5]*1e6;
[poincareSection.numPoints] = deal(100);
rOrbit = hypot(diff(poincareSection.endPosition(:,1)),diff(poincareSection.endPosition(:,2)));
poincareSection.orbitMaxLength = 2*(2*pi*rOrbit);
dThresh = 1e-3;

% Strainlines
strainlineMaxLength = 1e8;
strainlineLocalMaxDistance = 4*deltaX;
strainlineOdeSolverOptions = odeset('relTol',1e-4);

% Stretchlines
stretchlineMaxLength = 1e8;
stretchlineLocalMaxDistance = 8*deltaX;
stretchlineOdeSolverOptions = odeset('relTol',1e-4);

% Graphics properties
strainlineColor = 'r';
stretchlineColor = 'b';
lambdaLineColor = [0,.6,0];

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

%% Lambda-line LCSs
[etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'dThresh',dThresh);

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
strainlineLcs = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineOdeSolverOptions);

% Remove strainlines inside elliptic regions
strainlineLcs = remove_strain_in_elliptic(strainlineLcs,closedLambdaLine{1}{1}{end});
strainlineLcs = remove_strain_in_elliptic(strainlineLcs,closedLambdaLine{1}{2}{end});

% Plot hyperbolic strainline LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs,'UniformOutput',false);
hStrainlineLcs = [hStrainlineLcs{:}];
set(hStrainlineLcs,'color',strainlineColor)

uistack(hLambdaLineLcs,'top')
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(domain);
title(hAxes,'Stretchline and \lambda-line LCSs')

% Plot lambda-line LCSs
hLambdaLineLcs = copyobj(hLambdaLineLcs,hAxes);
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
stretchlineLcs = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions);

% Remove stretchlines inside elliptic regions
stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,closedLambdaLine{1}{1}{end});
stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,closedLambdaLine{1}{2}{end});

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs,'UniformOutput',false);
hStretchlineLcs = [hStretchlineLcs{:}];
set(hStretchlineLcs,'color',stretchlineColor)

uistack(hLambdaLineLcs,'top')
