%% Load ocean data set
load('ocean_geostrophic_velocity.mat');

% Set velocity to zero at boundaries
vLon(:,[1,end],:) = 0;
vLon(:,:,[1,end]) = 0;
vLat(:,[1,end],:) = 0;
vLat(:,:,[1,end]) = 0;

%% Set parameters
% Define right hand side of ODE, ocean.flow.derivative
interpMethod = 'spline';
vLon_interpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLat_interpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);
ocean.flow.derivative = @(t,y,~)derivative(t,y,vLon_interpolant,vLat_interpolant);

% Set domain of initial conditions
% Center of domain [lon,lat]
center = [3,-31];
halfwidth = 3;
ocean.flow.domain = [center(1)-halfwidth,center(1)+halfwidth;center(2)-halfwidth,center(2)+halfwidth];

% Use vectorized integration. 
ocean.flow.coupledIntegration = true;
% Set, if periodic boundary conditions in x and y direction
ocean.flow.periodicBc = [false,false];
% Set if incompressibility of the flow is enforced
ocean.flow.imposeIncompressibility = true;
% Set resolution of subdomain
nxy = 50;
ocean.flow.resolution = [nxy,nxy];

% Set integration time span (days)
ocean.flow.timespan = [98,128];

%% Animate ocean flow
ocean.flow = animate_flow(ocean.flow);

