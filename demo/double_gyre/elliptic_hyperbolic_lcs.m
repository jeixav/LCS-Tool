%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
timespan = [0,10];
domain = [0,2;0,1];
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% Lambda lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [.55,.55;.1,.1];
poincareSection(2).endPosition = [1.53,.45;1.95,.05];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end
lambda = .93:.01:1.07;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
forceEtaComplexNaN = true;

% Shrink lines
shrinkLineMaxLength = 20;
shrinkLineLocalMaxDistance = 2*deltaX;
shrinkLineOdeSolverOptions = odeset('relTol',1e-6);

% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 10*deltaX;
stretchLineOdeSolverOptions = odeset('relTol',1e-6);

% Graphic properties
repellingColor = 'r';
attractingColor = 'b';
ellipticColor = [0,.6,0];

hAxes = setup_figure(domain);
title(hAxes,'Repelling and elliptic LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

%% Elliptic LCSs
[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'odeSolverOptions',lambdaLineOdeSolverOptions);

ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)];

% Plot elliptic LCSs
hEllipticLcs = plot_elliptic_lcs(hAxes,ellipticLcs);
set(hEllipticLcs,'color',ellipticColor)
set(hEllipticLcs,'linewidth',2)
drawnow

%% Hyperbolic repelling LCSs
shrinkLine = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions);

% Remove shrink lines inside elliptic LCSs
for i = 1:nPoincareSection
    shrinkLine = remove_strain_in_elliptic(shrinkLine,ellipticLcs{i});
end

% Plot repelling LCSs
hRepellingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),shrinkLine,'UniformOutput',false);
hRepellingLcs = [hRepellingLcs{:}];
set(hRepellingLcs,'color',repellingColor)

uistack(hEllipticLcs,'top')
drawnow

%% Hyperbolic attracting LCSs
hAxes = setup_figure(domain);
title(hAxes,'Attracting and elliptic LCSs')

% Copy objects from repelling LCS plot
hEllipticLcs = copyobj(hEllipticLcs,hAxes);
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
stretchLine = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
end

% Plot attracting LCSs
hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchLine,'UniformOutput',false);
hAttractingLcs = [hAttractingLcs{:}];
set(hAttractingLcs,'color',attractingColor)

uistack(hEllipticLcs,'top')
