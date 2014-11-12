%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [0,lengthX;[-1,1]*2.25*lengthY];
timespan = [0,2*lengthX/u];
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
periodicBc = [true,false];

%% LCS parameters

% Lambda-lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.25e6,1.5e6;1.4e6,2.6e6];
poincareSection(2).endPosition = [6.5e6,-1.4e6;5e6,-3e6];
poincareSection(3).endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection(4).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection(5).endPosition = [1.65e7, 1.5e6;1.5e7, 2.6e6];
[poincareSection.numPoints] = deal(80);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end
lambda = .9:.01:1.1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
forceEtaComplexNaN = true;

% Strainlines
strainlineMaxLength = 1e8;
strainlineLocalMaxDistance = 8*deltaX;
strainlineOdeSolverOptions = odeset('relTol',1e-4);

% Stretchlines
stretchlineMaxLength = 1e8;
stretchlineLocalMaxDistance = 4*deltaX;
stretchlineOdeSolverOptions = odeset('relTol',1e-4);

% Graphics properties
repellingColor = 'r';
attractingColor = 'b';
ellipticColor = [0,.6,0];
lcsInitialPositionMarkerSize = 2;

hAxes = setup_figure(domain);
title(hAxes,'Repelling and elliptic LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);

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

[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'lambdaLineOdeSolverOptions',lambdaLineOdeSolverOptions,'periodicBc',periodicBc);

ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)];

% Plot elliptic LCSs
hEllipticLcs = plot_elliptic_lcs(hAxes,ellipticLcs);
set(hEllipticLcs,'color',ellipticColor)
set(hEllipticLcs,'linewidth',2)

% Plot closed lambda lines
hClosedLambdaLinePos = plot_closed_orbit(hAxes,closedLambdaLinePos);
hClosedLambdaLineNeg = plot_closed_orbit(hAxes,closedLambdaLineNeg);
hClosedLambdaLine = [hClosedLambdaLinePos,hClosedLambdaLineNeg];
set(hClosedLambdaLine,'color',ellipticColor)
drawnow

%% Hyperbolic repelling LCSs
[strainlineLcs,strainlineLcsInitialPosition] = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineOdeSolverOptions,'periodicBc',periodicBc);

% Remove strainlines inside elliptic LCSs
for i = 1:nPoincareSection
    strainlineLcs = remove_strain_in_elliptic(strainlineLcs,ellipticLcs{i});
    idx = inpolygon(strainlineLcsInitialPosition(1,:),strainlineLcsInitialPosition(2,:),ellipticLcs{i}(:,1),ellipticLcs{i}(:,2));
    strainlineLcsInitialPosition = strainlineLcsInitialPosition(:,~idx);
end

% Plot hyperbolic repelling LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs,'UniformOutput',false);
hStrainlineLcs = [hStrainlineLcs{:}];
set(hStrainlineLcs,'color',repellingColor)
hStrainlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineLcsInitialPosition(1,idx),strainlineLcsInitialPosition(2,idx)),1:size(strainlineLcsInitialPosition,2),'UniformOutput',false);
hStrainlineLcsInitialPosition = [hStrainlineLcsInitialPosition{:}];
set(hStrainlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStrainlineLcsInitialPosition,'marker','o')
set(hStrainlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineLcsInitialPosition,'MarkerFaceColor',repellingColor)

uistack(hEllipticLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
drawnow

%% Hyperbolic attracting LCSs
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
% minimums
[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions,'periodicBc',periodicBc);

% Remove stretchlines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,ellipticLcs{i});
    idx = inpolygon(stretchlineLcsInitialPosition(1,:),stretchlineLcsInitialPosition(2,:),ellipticLcs{i}(:,1),ellipticLcs{i}(:,2));
    stretchlineLcsInitialPosition = stretchlineLcsInitialPosition(:,~idx);
end

% Plot hyperbolic attracting LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs,'UniformOutput',false);
hStretchlineLcs = [hStretchlineLcs{:}];
set(hStretchlineLcs,'color',attractingColor)
hStretchlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineLcsInitialPosition(1,idx),stretchlineLcsInitialPosition(2,idx)),1:size(stretchlineLcsInitialPosition,2),'UniformOutput',false);
hStretchlineLcsInitialPosition = [hStretchlineLcsInitialPosition{:}];
set(hStretchlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStretchlineLcsInitialPosition,'marker','o')
set(hStretchlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchlineLcsInitialPosition,'MarkerFaceColor',attractingColor)

uistack(hEllipticLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
