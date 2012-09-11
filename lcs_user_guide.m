%% Lagrangian Coherent Structures Toolbox User Guide
%

%% Introduction
% The LCS toolbox is demonstrated by analyzing a double gyre flow:
% 
% $\dot x = -\pi A \sin[-\pi f(x,t)] \cos(\pi y)$
%
% $\dot y = \pi A \cos[\pi f(x,t)] \sin(\pi y) \frac{\partial f(x,t)}{\partial x}$
%
% where
%
% $f(x,t) = \epsilon \sin(\omega t)x^2 + [1 - 2 \epsilon \sin(\omega t)] x$
%
% and $A$, $\epsilon$ and $\omega$ are constants. Further details are
% availaibe in <http://dx.doi.org/10.1016/j.physd.2005.10.007 DOI:10.1016/j.physd.2005.10.007>, 
% <http://dx.doi.org/10.5194/npg-7-59-2000 DOI:10.5194/npg-7-59-2000>,
% and <http://dx.doi.org/10.5194/npg-4-223-1997 DOI:10.5194/npg-4-223-1997>.

%% Flow definition
% The flow vector field is defined as a symbolic function. The definition
% for the double gyre is:
t = sym('t');
x = sym('x');
y = sym('y');

p = struct('epsilon',.1,'a',.1,'omega',pi/5);

forcing = p.epsilon*sin(p.omega*t)*x^2 + (1 - 2*p.epsilon...
    *sin(p.omega*t))*x;

flow.symDerivative(1) = -pi*p.a*sin(pi*forcing)*cos(pi*y);
flow.symDerivative(2) = pi*p.a*cos(pi*forcing).*sin(pi*y)...
    *(2*p.epsilon*sin(p.omega*t)*x + 1 - 2*p.epsilon*sin(p.omega*t));

%%
% The flow domain, timespan and resolution must be defined also:
flow = set_flow_domain([0 2; 0 1],flow);
flow = set_flow_timespan([0 20],flow);
flow = set_flow_resolution(uint64([2 1]*10),flow);

%% Flow animation
% To verify that the flow has been correctly defined, it can be animated:
flow = animate_flow(flow);

%%
% Parameters can be changed and the animation re-run. For example
flow = set_flow_timespan([0 30],flow);
flow = set_flow_resolution(uint64([2 1]*20),flow);
flow = animate_flow(flow);

%% Hyperbolic barriers
% Hyperbolic barriers are obtained from strainlines. Strainlines are
% computed based on a resolution representing a grid of initial conditions:
strainline = set_strainline_resolution(uint64([2 1]*5));

%%
% A maximum length for strainlines must be specified. Strainlines are
% integrated until reaching the boundary. Nonetheless, a maximum length
% needs to be specified to bound integration time. This maximum length
% is found heuristically.
strainline = set_strainline_max_length(10,strainline);

%%
% The following parameters are used to filter LCSs from all calculated
% strainlines. To start, these parameters are set to display all
% strainlines.
strainline = set_strainline_geodesic_deviation_tol(inf,strainline);
strainline = set_strainline_length_tol(0,strainline);
strainline.filteringMethod = 'hausdorff';
strainline.filteringDistanceTol = 0;

%%
% This specifies everything necessary. The function strain_lcs_script is
% calculates and plots all strainlines:
doubleGyre = struct('flow',flow,'strainline',strainline);
doubleGyre = strain_lcs_script(doubleGyre);

%%
% The strainlines appear quite jagged. To fix this, the flow resolution
% needs to be increased and the ODE integration error tolerance decreased.
doubleGyre.flow = set_flow_resolution(uint64([2 1]*1000),doubleGyre.flow);
doubleGyre.strainline = set_strainline_ode_solver_options(...
    odeset('relTol',1e-6),doubleGyre.strainline);
load('datasets/doubleGyre') % Load precomputed data for speed
doubleGyre = strain_lcs_script(doubleGyre);

%%
% Filtering parameters are adjusted to find significant hyperbolic barriers
doubleGyre.strainline = set_strainline_geodesic_deviation_tol(.05,...
    doubleGyre.strainline);
doubleGyre.strainline = set_strainline_length_tol(.5,...
    doubleGyre.strainline);
doubleGyre.strainline = set_strainline_hausdorff_distance(.5,...
    doubleGyre.strainline);
doubleGyre = strain_lcs_script(doubleGyre);

%% Parabolic barriers
% Parabolic barriers are obtained from shearlines. The calculation of
% shearlines is similar to strainlines. A grid of initial conditions is
% specified:
doubleGyre.flow = set_flow_resolution(uint64([2 1]*50),doubleGyre.flow);
doubleGyre.shearline = set_shearline_resolution(uint64([2 1]*5));

%%
% A maximum length for shearlines is specified. This length is found
% heuristically
doubleGyre.shearline = set_shearline_max_length(10,doubleGyre.shearline);

%%
% The geodesic deviation is used as a filtering parameter
doubleGyre.shearline = set_shearline_average_geodesic_deviation_tol(...
    [inf inf],doubleGyre.shearline);

%%
% The ODE integration error tolerance should be decreased. Then the
% function shear_lcs_script calculates and plots all shearlines.
doubleGyre.shearline = set_strainline_ode_solver_options(...
    odeset('relTol',1e-6),doubleGyre.shearline);
doubleGyre = shear_lcs_script(doubleGyre);

%%
% Filtering parameters are adjusted to find significant shearlines
doubleGyre.shearline = set_shearline_average_geodesic_deviation_tol([1 1],doubleGyre.shearline);
doubleGyre = shear_lcs_script(doubleGyre);

%% Elliptic barriers
%
