%% Load ocean data set
velocityDataFile = fullfile('..','..','datasets','ocean_fhuhn','ocean_geostrophic_velocity.mat');
load(velocityDataFile);

% Set velocity to zero at boundaries
vlon(:,[1,end],:) = 0;
vlon(:,:,[1,end]) = 0;
vlat(:,[1,end],:) = 0;
vlat(:,:,[1,end]) = 0;

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

% Use vectorized integration. 
ocean.flow.coupledIntegration = true;
% Set, if periodic boundary conditions in x and y direction
ocean.flow.periodicBc = [false,false];
% Set if incompressibility of the flow is enforced
ocean.flow.imposeIncompressibility = true;
% Set resolution of subdomain
nxy = 50;
subdomainResolution = [nxy,nxy];
ocean.flow = set_flow_resolution(subdomainResolution,ocean.flow);

% Set integration time span (days)
ocean.flow = set_flow_timespan([98,128],ocean.flow);

%% Animate ocean flow
ocean.flow = animate_flow(ocean.flow);

