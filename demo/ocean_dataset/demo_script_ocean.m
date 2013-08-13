%% Demo script of LCS Tool
% Hyperbolic and elliptic LCS in an ocean data set
% fhuhn - 2013/07/25

close all
clear all
clc

addpath('../../');

%% Load ocean velocity data set
fprintf('Load velocity data ...\n');
velocityDataFile = '../../datasets/ocean_fhuhn/ocean_geostrophic_velocity.mat';
load(velocityDataFile);

% Define velocity domain
velDomain = [lon(1) lon(end); lat(1) lat(end)];

% Set velocity to zero at boundaries
vlon(:,1, 1:end) = 0;
vlon(:,end, 1:end) = 0;
vlon(:,1:end, 1) = 0;
vlon(:,1:end, end) = 0;
vlat(:,1, 1:end) = 0;
vlat(:,end, 1:end) = 0;
vlat(:,1:end, 1) = 0;
vlat(:,1:end, end) = 0;

% Control plot of instantaneous geostrophic velocity field
hFigVelData = figure; 
hAxesVelData = axes('parent',hFigVelData);
hold(hAxesVelData,'on');
hVelAbs = imagesc(lon, lat, sqrt(squeeze(vlon(1,:,:)).^2 + squeeze(vlat(1,:,:)).^2),'Parent',hAxesVelData);
hVelField = quiver(hAxesVelData,lon, lat, squeeze(vlon(1,:,:)), squeeze(vlat(1,:,:)), 2, 'k');
xlabel(hAxesVelData,'Lon [^\circ]'); ylabel(hAxesVelData,'Lat [^\circ]');
title(hAxesVelData,'Instantaneous geostrophic velocity field - |v| [degree/day]');
colorbar('peer',hAxesVelData);
axis(hAxesVelData,'equal','tight');
hold(hAxesVelData,'off');
drawnow

%% Set parameters
% Define right hand side of ODE, ocean.flow.derivative
interpMethod = 'spline';
vlon_interpolant = griddedInterpolant( {time, lat, lon}, vlon, interpMethod);
vlat_interpolant = griddedInterpolant( {time, lat, lon}, vlat, interpMethod);
ocean.flow.derivative = @(t,y,useEoV)flowdata_derivative(t,y,useEoV,vlon_interpolant,vlat_interpolant);

% Set domain of initial conditions
%*****************************
% Center of domain [lon lat]
center = [3.0 -31.0];
halfwidth = 3.0;
%*****************************
subdomain = [center(1)-halfwidth center(1)+halfwidth; center(2)-halfwidth center(2)+halfwidth]; % lon lat
ocean.flow = set_flow_domain(subdomain, ocean.flow);

% Set blocksize for vectorized integration. 
ocean.flow.coupledIntegration = true;
% Set, if periodic boundary conditions in x and y direction
ocean.flow.periodicBc = [false false];
% Set computation method for Cauchy-Green (CG) tensor
ocean.flow.cgStrainMethod.name = 'finiteDifference';
% Set if CG eigenvalues computed from main grid ('true' results in smoother eigenvalue fields)
ocean.flow.cgStrainMethod.eigenvalueFromMainGrid = false;
% Set auxiliary grid distance (relative value, i.e. 0.1 means 10% of maingrid size)
ocean.flow.cgStrainMethod.auxiliaryGridRelativeDelta = 0.1;
% Set computation method for eigenvectors
% false: use 'eig' function of MATLAB
% true: xi2 explicitly from auxiliary grid CG, xi1 as rotated xi2
ocean.flow.customEigMethod = true;
% Set if incompressibility of the flow is enforced,
%i.e., lambda1 = 1/lamda2
ocean.flow.imposeIncompressibility = true;
% Set resolution of subdomain
nxy = 400;
subdomainResolution = [nxy nxy];
ocean.flow = set_flow_resolution(subdomainResolution,ocean.flow);
% Define grid vectors of subdomain
lonAxis = linspace(ocean.flow.domain(1,1), ocean.flow.domain(1,2), ocean.flow.resolution(1));
latAxis = linspace(ocean.flow.domain(2,1), ocean.flow.domain(2,2), ocean.flow.resolution(2));

gridSpace = diff(ocean.flow.domain(1,:))/(double(ocean.flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;
ocean.strainline = set_strainline_max_length(20);
ocean.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),ocean.strainline);

%% **********************************
% Repelling LCS - forward in time
% *********************************** 

% Set integration time span (days)
ocean.flow = set_flow_timespan([98 128],ocean.flow);

%% Compute Cauchy-Green strain eigenvalues and eigenvectors
fprintf('Integrate flow forward ...\n');
[ocean.flow.cgEigenvalue,ocean.flow.cgEigenvector] = eig_cgStrain(ocean.flow,ocean.flow.cgStrainMethod,ocean.flow.customEigMethod,ocean.flow.coupledIntegration);
cgEigenvalue = reshape(ocean.flow.cgEigenvalue,[fliplr(ocean.flow.resolution),2]);
cgEigenvector = reshape(ocean.flow.cgEigenvector,[fliplr(ocean.flow.resolution),4]);

%% SHEARLINES
[ocean.shearline.etaPos,ocean.shearline.etaNeg] = lagrangian_shear(ocean.flow.cgEigenvector,ocean.flow.cgEigenvalue);

%% Define poincare sections for closed orbit detection
% Poincare section should be placed with 1st point in center of elliptic region and
% with second point outside the elliptic region
%
% A poincare section extending far into the hyperbolic region might cause
% long runtime for the shearline integration
%
% poincareSection{i}.endPosition = [x1 y1; x2 y2]
% poincareSection{i}.numPoints = nPoints
% poincareSection{i}.integrationLength = integrationLength

clear poincareSection

% Poincare section 1
poincareSection{1}.endPosition = [3.15 -32.2; 3.7 -31.6];
% Poincare section 2
poincareSection{2}.endPosition = [5.0 -31.6;  5.3 -31.6];
% Poincare section 3
poincareSection{3}.endPosition = [4.8 -29.5;  4.4 -29.5];
% Poincare section 4
poincareSection{4}.endPosition = [1.5 -30.9;  1.9 -31.1];
% Poincare section 5
poincareSection{5}.endPosition = [2.9 -29.2;  3.2 -29.0];
% Number of points along each poincare section
numPoints = 100;

nPoincareSection = numel(poincareSection);
for i=1:nPoincareSection
    poincareSection{i}.numPoints = numPoints;   %#ok<SAGROW>
    % radius = length of poincare section
    rOrbit = norm(poincareSection{i}.endPosition(2,:)-poincareSection{i}.endPosition(1,:));
    % Set integration length conservatively = twice the expected circumference
    poincareSection{i}.integrationLength = [0 2*(2*pi*rOrbit)]; %#ok<SAGROW>
end

%% Closed orbit detection
%*************************
odeSolverOptions = odeset('relTol',1e-6);
% Set to 'true' to see return maps of poincare sections
showPlot = false;
dThresh = 1e-2;
%*************************
fprintf('Detect elliptic LCS ...\n');
[closedOrbits, orbits] = poincare_closed_orbit_multi(ocean.flow,ocean.shearline,poincareSection,odeSolverOptions,dThresh,showPlot); %#ok<NASGU>

%% Plot elliptic LCS
% Set up figure
hAxes = setup_figure(ocean.flow.domain);
colormap(hAxes,'gray');
xlabel(hAxes,'Lon [^\circ]'); ylabel(hAxes,'Lat [^\circ]');
% etaPos closed orbits
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}(:,1),closedOrbits{i}{1}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
% etaNeg closed orbits
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}(:,1),closedOrbits{i}{2}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
drawnow

%% Compute strainlines
fprintf('Detect hyperbolic LCS ...\n');
fprintf('Compute strainlines ...\n');
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,ocean.strainline.maxLength,...
    cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),ocean.flow.domain,ocean.flow.periodicBc);

for i=1:nPoincareSection
% Remove strainlines inside of ellitpic regions
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1});
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2});
% Set initial conditions to NaN inside of elliptic regions
idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{1}(:,1),closedOrbits{i}{1}(:,2));
strainlineInitialPosition(1,idx) = NaN;
strainlineInitialPosition(2,idx) = NaN;
idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{2}(:,1),closedOrbits{i}{2}(:,2));
strainlineInitialPosition(1,idx) = NaN;
strainlineInitialPosition(2,idx) = NaN;
end

%% Plot hyperbolic LCS
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','r');
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top'); 

%% Plot additional technical details of LCS detection
% To plot, set respective visible flag to 'on'

print_epspdfpng(gcf,'ocean_LCS_fwd');

% Plot FTLE in the background
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(ocean.flow.timespan));
hFtle = plot_ftle(hAxes,ocean.flow,ftle);
colormap(hAxes,flipud(gray));
uistack(hClosedOrbitsEtaNeg,'bottom');
uistack(hStrainline,'top');
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top');
set(hFtle,'visible','on');
% colorbar('hide');

% Plot poincare sections
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1), poincareSection{idx}.endPosition(:,2), 'ro-', 'linewidth',2), 1:size(poincareSection,2) );
set(hPoincareSection,'color','w');
set(hPoincareSection,'marker','o');
set(hPoincareSection,'markerFaceColor','w');
set(hPoincareSection,'linewidth',2);
drawnow
set(hPoincareSection,'visible','on');

% Plot initial conditions of strainlines at local maxima of lambda2
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'MarkerSize',5)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor','r')
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top');
uistack(hPoincareSection,'top');
set(hStrainlineInitialPosition,'visible','on');

% Plot all orbits at poincare sections
hOrbits = hggroup;
for j=1:nPoincareSection
    % etaPos.
    hOrbitsPos = arrayfun(@(i)plot(hAxes,orbits{j}{1}{i}(:,1),orbits{j}{1}{i}(:,2),'b-'),1:numPoints);
    set(hOrbitsPos,'Parent',hOrbits);
    % etaNeg
    hOrbitsNeg = arrayfun(@(i)plot(hAxes,orbits{j}{2}{i}(:,1),orbits{j}{2}{i}(:,2),'k-'),1:numPoints);
    set(hOrbitsNeg,'Parent',hOrbits);
end
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top');
uistack(hPoincareSection,'top');
set(hOrbits,'visible','on');

print_epspdfpng(gcf,'ocean_LCS_fwd_details');

%% ********************************** 
% Attracting LCS - backward in time
% ***********************************

% Set integration time span (days)
ocean.flow = set_flow_timespan([98 68],ocean.flow);

%% Compute Cauchy-Green strain eigenvalues and eigenvectors
fprintf('Integrate flow backward ...\n');
[ocean.flow.cgEigenvalue,ocean.flow.cgEigenvector] = eig_cgStrain(ocean.flow,ocean.flow.cgStrainMethod,ocean.flow.customEigMethod,ocean.flow.coupledIntegration);
cgEigenvalue = reshape(ocean.flow.cgEigenvalue,[fliplr(ocean.flow.resolution),2]);
cgEigenvector = reshape(ocean.flow.cgEigenvector,[fliplr(ocean.flow.resolution),4]);

%% SHEARLINES
[ocean.shearline.etaPos,ocean.shearline.etaNeg] = lagrangian_shear(ocean.flow.cgEigenvector,ocean.flow.cgEigenvalue);

%% Define poincare sections for closed orbit detection
% Poincare section should be placed with 1st point in center of elliptic region and
% with second point outside the elliptic region
%
% A poincare section extending far into the hyperbolic region might cause
% long runtime for the shearline integration
%
% poincareSection{i}.endPosition = [x1 y1; x2 y2]
% poincareSection{i}.numPoints = nPoints
% poincareSection{i}.integrationLength = integrationLength

clear poincareSection

% Poincare section 1
poincareSection{1}.endPosition = [3.15 -32.2; 3.7 -31.6];
% Poincare section 2
poincareSection{2}.endPosition = [4.8 -29.5;  4.3 -29.7];
% Number of points along each poincare section
numPoints = 100;

nPoincareSection = numel(poincareSection);
for i=1:nPoincareSection
    poincareSection{i}.numPoints = numPoints;   %#ok<SAGROW>
    % radius = length of poincare section
    rOrbit = norm(poincareSection{i}.endPosition(2,:)-poincareSection{i}.endPosition(1,:));
    % Set integration length conservatively = twice the expected circumference
    poincareSection{i}.integrationLength = [0 2*(2*pi*rOrbit)]; %#ok<SAGROW>
end

%% Closed orbit detection
%*************************
odeSolverOptions = odeset('relTol',1e-6);
% Set to 'true' to see return maps of poincare sections
showPlot = false;
dThresh = 1e-2;
%*************************
fprintf('Detect elliptic LCS ...\n');
[closedOrbits, orbits] = poincare_closed_orbit_multi(ocean.flow,ocean.shearline,poincareSection,odeSolverOptions,dThresh,showPlot);

%% Plot elliptic LCS
% Set up figure
hAxes = setup_figure(ocean.flow.domain);
colormap(hAxes,'gray');
xlabel(hAxes,'Lon [^\circ]'); ylabel(hAxes,'Lat [^\circ]');
% etaPos closed orbits
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}(:,1),closedOrbits{i}{1}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
% etaNeg closed orbits
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}(:,1),closedOrbits{i}{2}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
drawnow

%% Compute strainlines
fprintf('Detect hyperbolic LCS ...\n');
fprintf('Compute strainlines ...\n');
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,ocean.strainline.maxLength,...
    cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),ocean.flow.domain,ocean.flow.periodicBc);

for i=1:nPoincareSection
% Remove strainlines inside of ellitpic regions
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1});
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2});
% Set initial conditions to NaN inside of elliptic regions
idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{1}(:,1),closedOrbits{i}{1}(:,2));
strainlineInitialPosition(1,idx) = NaN;
strainlineInitialPosition(2,idx) = NaN;
idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{2}(:,1),closedOrbits{i}{2}(:,2));
strainlineInitialPosition(1,idx) = NaN;
strainlineInitialPosition(2,idx) = NaN;
end

%% Plot hyperbolic LCS
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','b');
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top'); 

%% Plot additional technical details of LCS detection
% To plot, set respective visible flag to 'on'

print_epspdfpng(gcf,'ocean_LCS_bwd');

% Plot FTLE in the background
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(ocean.flow.timespan));
hFtle = plot_ftle(hAxes,ocean.flow,ftle);
uistack(hClosedOrbitsEtaNeg,'bottom');
uistack(hStrainline,'top');
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top');
set(hFtle,'visible','on'); 
% colorbar('hide');

% Plot poincare sections
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1), poincareSection{idx}.endPosition(:,2), 'ro-', 'linewidth',2), 1:size(poincareSection,2) );
set(hPoincareSection,'color','w');
set(hPoincareSection,'marker','o');
set(hPoincareSection,'markerFaceColor','w');
set(hPoincareSection,'linewidth',2);
set(hPoincareSection,'visible','on');

% Plot initial conditions of strainlines at local maxima of lambda2
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'MarkerSize',5)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor','b')
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top');
uistack(hPoincareSection,'top');
set(hStrainlineInitialPosition,'visible','on');

% Plot all orbits at poincare sections
hOrbits = hggroup;
for j=1:nPoincareSection
    % etaPos.
    hOrbitsPos = arrayfun(@(i)plot(hAxes,orbits{j}{1}{i}(:,1),orbits{j}{1}{i}(:,2),'r-'),1:numPoints);
    set(hOrbitsPos,'Parent',hOrbits);
    % etaNeg
    hOrbitsNeg = arrayfun(@(i)plot(hAxes,orbits{j}{2}{i}(:,1),orbits{j}{2}{i}(:,2),'k-'),1:numPoints);
    set(hOrbitsNeg,'Parent',hOrbits);
end
uistack(hClosedOrbitsEtaPos,'top');
uistack(hClosedOrbitsEtaNeg,'top');
uistack(hPoincareSection,'top');
set(hOrbits,'visible','on');

print_epspdfpng(gcf,'ocean_LCS_bwd_details');

%%


