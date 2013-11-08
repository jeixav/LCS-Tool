%% Load data set
load('ocean_geostrophic_velocity.mat')

% Set velocity to zero at boundaries
vlon(:,[1,end],:) = 0;
vlon(:,:,[1,end]) = 0;
vlat(:,[1,end],:) = 0;
vlat(:,:,[1,end]) = 0;

forwardLcsColor = 'r';
backwardLcsColor = 'b';
shearLcsColor = [0,.6,0];

%% Set parameters
% Define right hand side of ODE, ocean.flow.derivative
interpMethod = 'spline';
vlon_interpolant = griddedInterpolant({time,lat,lon},vlon,interpMethod);
vlat_interpolant = griddedInterpolant({time,lat,lon},vlat,interpMethod);
ocean.flow.derivative = @(t,y,~)flowdata_derivative(t,y,vlon_interpolant,vlat_interpolant);

% Set domain of initial conditions
% Center of domain [lon,lat]
center = [3,-31];
halfwidth = 3;
subdomain = [center(1)-halfwidth,center(1)+halfwidth;center(2)-halfwidth,center(2)+halfwidth];
ocean.flow = set_flow_domain(subdomain,ocean.flow);

% Set computation method for Cauchy-Green (CG) tensor
ocean.flow.cgStrainMethod.name = 'finiteDifference';
% Set if CG eigenvalues computed from main grid ('true' yields smoother eigenvalue fields)
ocean.flow.cgStrainMethod.eigenvalueFromMainGrid = false;
% Set auxiliary grid distance (relative value, i.e. 0.1 means 10% of maingrid size)
ocean.flow.cgStrainMethod.auxiliaryGridRelativeDelta = 0.1;
% Set computation method for eigenvectors
% false: use 'eig' function of MATLAB
% true: xi2 explicitly from auxiliary grid CG, xi1 as rotated xi2
ocean.flow.customEigMethod = false;
% Set if incompressibility of the flow is enforced,
% i.e., lambda1 = 1/lamda2
ocean.flow.imposeIncompressibility = true;
% Set resolution of subdomain
nxy = 400;
subdomainResolution = [nxy,nxy];
ocean.flow = set_flow_resolution(subdomainResolution,ocean.flow);

shearlineOdeSolverOptions = odeset('relTol',1e-6);
strainlineOdeSolverOptions = odeset('relTol',1e-4);

gridSpace = diff(ocean.flow.domain(1,:))/(double(ocean.flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;
ocean.strainline = set_strainline_max_length(20);

%% Repelling LCS - forward in time

% Set integration time span (days)
ocean.flow = set_flow_timespan([98,128],ocean.flow);

% Compute Cauchy-Green strain eigenvalues and eigenvectors
disp('Integrate flow forward ...')
[ocean.flow.cgEigenvalue,ocean.flow.cgEigenvector] = eig_cgStrain(ocean.flow,ocean.flow.cgStrainMethod,ocean.flow.customEigMethod);
cgEigenvalue = reshape(ocean.flow.cgEigenvalue,[fliplr(ocean.flow.resolution),2]);
cgEigenvector = reshape(ocean.flow.cgEigenvector,[fliplr(ocean.flow.resolution),4]);

% Shearlines
lambda = 1;
[ocean.shearline.etaPos,ocean.shearline.etaNeg] = lambda_line(ocean.flow.cgEigenvector,ocean.flow.cgEigenvalue,lambda);

% Define Poincare sections for closed orbit detection
% Poincare section should be placed with 1st point in center of elliptic region and
% with second point outside the elliptic region

poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

% poincareSection(i).endPosition = [longitude1,latitude1;longitude2,latitude2]
poincareSection(1).endPosition = [3.15,-32.2;3.7,-31.6];
poincareSection(2).endPosition = [5,-31.6;5.3,-31.6];
poincareSection(3).endPosition = [4.8,-29.5;4.4,-29.5];
poincareSection(4).endPosition = [1.5,-30.9;1.9,-31.1];
poincareSection(5).endPosition = [2.9,-29.2;3.2,-29];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length conservatively to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

% Closed orbit detection
disp('Detect elliptic LCS ...')
closedOrbits = poincare_closed_orbit_multi(ocean.flow,ocean.shearline,poincareSection,'odeSolverOptions',shearlineOdeSolverOptions);

% Plot elliptic LCS
hAxes = setup_figure(ocean.flow.domain);
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')
title(hAxes,'Forward-time LCS')
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbit
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)
drawnow

% Compute strainlines
disp('Detect hyperbolic LCS ...')
disp('Compute strainlines ...')
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,ocean.strainline.maxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),ocean.flow.domain,'odeSolverOptions',strainlineOdeSolverOptions);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
end

% Plot hyperbolic LCS
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',forwardLcsColor)
uistack(hClosedOrbitsEtaPos,'top')
uistack(hClosedOrbitsEtaNeg,'top')
drawnow

%% Attracting LCS - backward in time

% Set integration time span (days)
ocean.flow = set_flow_timespan([98,68],ocean.flow);

% Compute Cauchy-Green strain eigenvalues and eigenvectors
disp('Integrate flow backward ...')
[ocean.flow.cgEigenvalue,ocean.flow.cgEigenvector] = eig_cgStrain(ocean.flow,ocean.flow.cgStrainMethod,ocean.flow.customEigMethod);
cgEigenvalue = reshape(ocean.flow.cgEigenvalue,[fliplr(ocean.flow.resolution),2]);
cgEigenvector = reshape(ocean.flow.cgEigenvector,[fliplr(ocean.flow.resolution),4]);

% Shearlines
[ocean.shearline.etaPos,ocean.shearline.etaNeg] = lambda_line(ocean.flow.cgEigenvector,ocean.flow.cgEigenvalue,lambda);

% Poincare sections
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6];
poincareSection(2).endPosition = [4.8,-29.5;4.3,-29.7];

% Number of points along each Poincare section
[poincareSection.numPoints] = deal(100);

nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

% Closed orbit detection
shearlineOdeSolverOptions = odeset('relTol',1e-6);
showPoincareReturnMap = false;
dThresh = 1e-2;
disp('Detect elliptic LCS ...')
closedOrbits = poincare_closed_orbit_multi(ocean.flow,ocean.shearline,poincareSection,'odeSolverOptions',shearlineOdeSolverOptions);

% Plot elliptic LCS
hAxes = setup_figure(ocean.flow.domain);
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')
title(hAxes,'Backward-time LCS')
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbits
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)
drawnow

% Compute strainlines
disp('Detect hyperbolic LCS ...')
disp('Compute strainlines ...')
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,ocean.strainline.maxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),ocean.flow.domain,'odeSolverOptions',strainlineOdeSolverOptions);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
end

% Plot hyperbolic LCS
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',backwardLcsColor)
uistack(hClosedOrbitsEtaPos,'top')
uistack(hClosedOrbitsEtaNeg,'top')
