%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];
domain = [2e6,.55*lengthX;[-1,.25]*2.25*lengthY];
timespan = [0,2*lengthX/u];

% Make x and y grid spacing as equal as possible
resolutionX = 500;
gridSpace = diff(domain(1,:))/(double(resolutionX)-1);
resolutionY = round(diff(domain(2,:))/gridSpace);
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
lambdaStep = 0.02;
lambdaRange = 0.90:lambdaStep:1.10;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
poincareSection.endPosition = [6.5,-1.4;4.5,-3.5]*1e6;
[poincareSection.numPoints] = deal(100);
rOrbit = hypot(diff(poincareSection.endPosition(:,1)),diff(poincareSection.endPosition(:,2)));
poincareSection.orbitMaxLength = 2*(2*pi*rOrbit);
dThresh = 1e-4;
nPoincareSection = numel(poincareSection);

% Strainlines
strainlineMaxLength = 1e8;
strainlineLocalMaxDistance = 4*gridSpace;
strainlineOdeSolverOptions = odeset('relTol',1e-4);

% Stretchlines
stretchlineMaxLength = 1e8;
stretchlineLocalMaxDistance = 8*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-4);

% Graphics properties
strainlineColor = 'r';
stretchlineColor = 'b';
lambdaLineColor = [0,.6,0];

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);

%% Lambda-line LCSs
% find closed orbits for range of lambda values
closedLambdaLineArea = zeros(1,nPoincareSection);
lambda0 = nan(1,nPoincareSection);
orbitArea = nan(1,2);
closedLambdaLine = cell(1,nPoincareSection);
k=0;
for lambda = lambdaRange
    k=k+1;
    
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
    shearline.etaPos = real(shearline.etaPos);
    shearline.etaNeg = real(shearline.etaNeg);      
    
    % find closed orbits
    closedLambdaLineCandidate = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions);
    
    % keep outermost closed orbit
    for i = 1:nPoincareSection
        for j=1:2 % etaPos,etaNeg
            orbitArea(j) = polyarea(closedLambdaLineCandidate{i}{j}{end}(:,1),closedLambdaLineCandidate{i}{j}{end}(:,2));
        end        
        if max(orbitArea) > closedLambdaLineArea(i)
            closedLambdaLineArea(i) = max(orbitArea);
            closedLambdaLine{i}{1}{1} = closedLambdaLineCandidate{i}{1}{1};
            closedLambdaLine{i}{2}{1} = closedLambdaLineCandidate{i}{2}{1};
            % keep lambda values associated to closed orbits
            lambda0(i) = lambda;
        end        
    end    
end

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
strainlineLcs = remove_strain_in_shear(strainlineLcs,closedLambdaLine{1}{1}{end});
strainlineLcs = remove_strain_in_shear(strainlineLcs,closedLambdaLine{1}{2}{end});

% Plot hyperbolic strainline LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs);
set(hStrainlineLcs,'color',strainlineColor)

uistack(hLambdaLineLcs,'top')
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(domain);
title(hAxes,'Stretchline and \lambda-line LCSs')

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

% Remove stretchlines inside elliptic regions
stretchlineLcs = remove_strain_in_shear(stretchlineLcs,closedLambdaLine{1}{1}{end});
stretchlineLcs = remove_strain_in_shear(stretchlineLcs,closedLambdaLine{1}{2}{end});

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs);
set(hStretchlineLcs,'color',stretchlineColor)

uistack(hLambdaLineLcs,'top')
