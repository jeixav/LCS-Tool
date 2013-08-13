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
% !!! USE A SQUARE DOMAIN !!!
%*****************************
center = [3.0 -31.0];
halfwidth = 3.0;
%*****************************
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
% FIXME Double gyre demos have a separate script only for flow animation.
% Consider doing the same with the ocean dataset, since flow animation
% requires the specification of a smaller number of parameters and a lower
% resolution.
% FIXME Should extend animate_flow function so xlabel and ylabel can be set
% to longitude and latitude
animationFilename = 'oceanFlowAnimation';
animate_flow(ocean.flow, animationTime, framerate, animationFilename);