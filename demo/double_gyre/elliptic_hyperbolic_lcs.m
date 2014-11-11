%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
timespan = [0,5];
resolutionX = 750;
resolutionY = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% Lambda-lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [.55,.55;.2,.5];
poincareSection(2).endPosition = [1.53,.45;1.9,.5];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);

% Strainlines
strainlineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLocalMaxDistance = 2*gridSpace;
strainlineOdeSolverOptions = odeset('relTol',1e-6);

% Stretchlines
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 10*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-6);

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
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions);

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
for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlineLcs = remove_strain_in_elliptic(strainlineLcs,closedLambdaLine{i}{1}{end});
    strainlineLcs = remove_strain_in_elliptic(strainlineLcs,closedLambdaLine{i}{2}{end});
end

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
for i = 1:nPoincareSection
    % Remove stretchlines inside elliptic regions
    stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,closedLambdaLine{i}{1}{end});
    stretchlineLcs = remove_strain_in_elliptic(stretchlineLcs,closedLambdaLine{i}{2}{end});
end

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs,'UniformOutput',false);
hStretchlineLcs = [hStretchlineLcs{:}];
set(hStretchlineLcs,'color',stretchlineColor)

uistack(hLambdaLineLcs,'top')
