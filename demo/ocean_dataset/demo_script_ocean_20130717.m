%% Demo script of LCS Tool
% Elliptic LCS in an ocean data set
% fhuhn - 2013/07/12

%% Set path variables
% FIXME Let user decide whether they want to close all their figures and
% clear all their variables. Do not force it inside this script.
close all
clear all
clc

% add LCS Tool folder
% FIXME Let user take care of adding path to LCS Tool once upon starting
% MATLAB, do not call addpath every time script is run.
addpath('../../');

% Open matlab pool for parallel computing
% FIXME Let user decide is they want to run with multiple workers. Do not
% open the matlabpool inside this demo.
if ~matlabpool('size')
    matlabpool('open');
end
% FIXME This command is no longer necessary, remove it
pctRunOnAll javaaddpath('../../ParforProgress2');

%% Load ocean velocity data set
% FIXME Please store datasets under ../../datasets/ocean_fhuhn
velocityDataFile = 'ocean_geostrophic_velocity.mat';
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
% FIXME Use graphics objects handles in commands below. For example
% hFigure = figure;
% hAxes = axes('parent',hFigure)
% hold(hAxes,on)
% ...
figure; hold on;
imagesc(lon, lat, sqrt(squeeze(vlon(1,:,:)).^2 + squeeze(vlat(1,:,:)).^2) );
quiver(lon, lat, squeeze(vlon(1,:,:)), squeeze(vlat(1,:,:)), 2, 'k')
xlabel('Lon [^\circ]'); ylabel('Lat [^\circ]');
title('Instantaneous geostrophic velocity field - |v| [degree/day]');
colorbar
axis equal tight
drawnow

%% Set parameters
% Define right hand side of ODE, ocean.flow.derivative
interpMethod = 'spline';
vlon_interpolant = griddedInterpolant( {time, lat, lon}, vlon, interpMethod);
vlat_interpolant = griddedInterpolant( {time, lat, lon}, vlat, interpMethod);
ocean.flow.derivative = @(t,y,useEoV)flowdata_derivative(t,y,useEoV,vlon_interpolant,vlat_interpolant);

% Set domain of initial conditions
% !!! USE A SQUARE DOMAIN !!!
%*****************************
center = [3.0 -31.0];
%*****************************
halfwidth = 3.0;
subdomain = [center(1)-halfwidth center(1)+halfwidth; center(2)-halfwidth center(2)+halfwidth]; % lon lat
ocean.flow = set_flow_domain(subdomain, ocean.flow);

% Set integration time span (days)
ocean.flow = set_flow_timespan([0 90],ocean.flow);

% Set blocksize for vectorized integration. 
ocean.flow.coupledIntegration = true;

% Set, if periodic boundary conditions in x and y direction
ocean.flow.periodicBc = [false false];

% Set computation method for Cauchy-Green (CG) tensor
ocean.flow.cgStrainMethod.name = 'finiteDifference';
% Set if CG eigenvalues computed from main grid (true results in smoother eigenvalue fields)
ocean.flow.cgStrainMethod.eigenvalueFromMainGrid = false;
% Set auxiliary grid distance (relative value, i.e. 0.1 means 10% of maingrid size)
ocean.flow.cgStrainMethod.auxiliaryGridRelativeDelta = 0.1;

% Set computation method for eigenvectors
ocean.flow.cgStrainEigMethod = 'custom';
% standard: use 'eig' function of MATLAB
% custom: xi2 explicitly from auxiliary grid CG, xi1 as rotated xi2

% Set if incompressibility of the flow is enforced
ocean.flow.imposeIncompressibility = true;
%i.e., lambda1 = 1/lamda2

%% Animate flow and save movie file
% set resolution of subdomain to low value for flow visualization
nxy = 50;
subdomainResolution = [nxy nxy];
ocean.flow = set_flow_resolution(subdomainResolution,ocean.flow);

% Set ode solver options for flow integration
odeOptions = odeset('RelTol',1e-6,'AbsTol',1e-8);
ocean.flow = set_flow_ode_solver_options(odeOptions,ocean.flow);

% Save .avi video file
framerate = 10;
animationTime = 10;
animationFilename = 'oceanFlowAnimation';
animate_flow(ocean.flow, animationTime, framerate, animationFilename);

%% Set resolution of subdomain to high resolution to detect LCS
nxy = 400;
subdomainResolution = [nxy nxy];
ocean.flow = set_flow_resolution(subdomainResolution,ocean.flow);
% Define grid vectors of subdomain
lonAxis = linspace(ocean.flow.domain(1,1), ocean.flow.domain(1,2), ocean.flow.resolution(1));
latAxis = linspace(ocean.flow.domain(2,1), ocean.flow.domain(2,2), ocean.flow.resolution(2));

%% STRAINLINES

% XXX

%% SHEARLINES
% set shearline parameters
ocean.shearline = set_shearline_resolution([1 1]);
ocean.shearline = set_shearline_max_length(10 ,ocean.shearline);
ocean.shearline = set_shearline_ode_solver_options(odeset('relTol',1e-6),ocean.shearline);
ocean.shearline = set_shearline_average_geodesic_deviation_tol([0.01 0.01],ocean.shearline);

ocean = shear_lcs_script(ocean);

%% save/load
% save run01_400x400_AuxGrid
% load run01_400x400_AuxGrid

%% Define CG eigenvalues/vectors
lambda1 = reshape(ocean.flow.cgEigenvalue(:,1), fliplr(ocean.flow.resolution));
lambda2 = reshape(ocean.flow.cgEigenvalue(:,2), fliplr(ocean.flow.resolution));
xi1x = reshape(ocean.flow.cgEigenvector(:,1), fliplr(ocean.flow.resolution));
xi1y = reshape(ocean.flow.cgEigenvector(:,2), fliplr(ocean.flow.resolution));
% xi2x = reshape(ocean.flow.cgEigenvector(:,3), fliplr(ocean.flow.resolution));
% xi2y = reshape(ocean.flow.cgEigenvector(:,4), fliplr(ocean.flow.resolution));

% Control plot of largest CG eigenvalue lambda2
figure;
imagesc(lonAxis, latAxis, log10(sqrt(lambda2)));
axis equal
axis tight
xlabel('Lon'); ylabel('Lat');
set(gca, 'ydir', 'normal');

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

% Poincare section 1
poincareSection{1}.endPosition = [3.15 -32.2; 3.7 -31.7];
% Poincare section 2
poincareSection{2}.endPosition = [5.0 -31.6;  5.4 -31.6];
% Poincare section 3
poincareSection{3}.endPosition = [1.5 -30.9;  2.0 -31.2];
% poincareSection{3}.endPosition = [5.2 -33.8;  4.6 -33.8];
% Poincare section 4
poincareSection{4}.endPosition = [4.8 -29.5;  4.4 -29.5];
% Number of points along each poincare section
numPoints = 100;

nPoincareSection = size(poincareSection,2);
for i=1:nPoincareSection
    poincareSection{i}.numPoints = numPoints;   %#ok<SAGROW>
    % radius = length of poincare section
    rOrbit = norm(poincareSection{i}.endPosition(2,:)-poincareSection{i}.endPosition(1,:));
    % Set integration length conservatively = twice the expected circumference
    poincareSection{i}.integrationLength = [0 2*(2*pi*rOrbit)]; %#ok<SAGROW>
end

figure; hold on;
imagesc(lonAxis, latAxis, log10(sqrt(lambda2)));
arrayfun(@(idx)plot(poincareSection{idx}.endPosition(:,1), poincareSection{idx}.endPosition(:,2), 'ro-', 'linewidth',2), 1:size(poincareSection,2) );
axis equal tight

%% Closed orbit detection
clc
close all

%*************************
odeSolverOptions = odeset('relTol',1e-6);
showPlot = 1;
dThresh = 1e-2;
%*************************
[closedOrbits, orbits] = poincare_closed_orbit_multi(ocean.flow,ocean.shearline,poincareSection,odeSolverOptions,dThresh,showPlot);

%% Plot results of closed orbit detection
hfigure = figure;
hAxes = axes;
hold on
imagesc(lonAxis,latAxis,log10(sqrt(lambda2)),'Parent',hAxes);
% Poincare sections
arrayfun(@(i)plot(hAxes,poincareSection{i}.endPosition(:,1),poincareSection{i}.endPosition(:,2),'ro--','linewidth',2),1:nPoincareSection);
% All orbits
for j=1:nPoincareSection
    % etaPos.
    arrayfun(@(i)plot(hAxes,orbits{j}{1}{i}(:,1),orbits{j}{1}{i}(:,2),'b-'),1:numPoints);
    % etaNeg
    arrayfun(@(i)plot(hAxes,orbits{j}{2}{i}(:,1),orbits{j}{2}{i}(:,2),'k-'),1:numPoints);
end
% etaPos closed orbits
arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}(:,1),closedOrbits{i}{1}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
% etaNeg closed orbits
arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}(:,1),closedOrbits{i}{2}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
colormap gray
axis equal tight
xlabel('Lon [^\circ]'); ylabel('Lat [^\circ]');

%%





