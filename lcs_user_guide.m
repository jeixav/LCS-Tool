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
% The ODE system for the double gyre is defined in a file called
% double_gyre_derivative.m. This file is in the flow_templates directory.
% To use this file, execute the following:
addpath flow_templates
flow.derivative = @(t,x)double_gyre_derivative(t,x);
flow.coupledIntegration = true;

%%
% The flow domain, timespan and resolution must be defined also:
flow = set_flow_domain([0 2; 0 1],flow);
flow = set_flow_timespan([0 20],flow);
flow = set_flow_resolution([2 1]*10,flow);

%% Flow animation
% To verify that the flow has been correctly defined, it can be animated:
flow = animate_flow(flow);

%%
% Parameters can be changed and the animation re-run. For example
flow = set_flow_timespan([0 30],flow);
flow = set_flow_resolution([2 1]*20,flow);
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
strainline = set_strainline_max_length(5,strainline);

%%
% The following parameters are used to filter LCSs from all calculated
% strainlines. To start, these parameters are set to display all
% strainlines.
strainline = set_strainline_geodesic_deviation_tol(inf,strainline);
strainline = set_strainline_length_tol(0,strainline);
strainline = set_strainline_filtering_method('superminimize',strainline);
filteringParameters.distance = 0;
filteringParameters.resolution = [1 1];
strainline = set_strainline_filtering_parameters(filteringParameters,strainline);

%%
% This specifies everything necessary. The function strain_lcs_script is
% calculates and plots all strainlines:
doubleGyre = struct('flow',flow,'strainline',strainline);
[doubleGyre,hAxes] = strain_lcs_script(doubleGyre);

%%
% The plot produced shows strainlines that meet filtering criteria. To see
% all the strainlines, execute:
set(findobj(hAxes,'tag','strainline'),'visible','on')

%%
% To get a list of all graphics objects whose visibility can be controlled,
% type:
unique(get(get(hAxes,'children'),'tag'))

%%
% The strainlines appear quite jagged. To fix this, the flow resolution
% needs to be increased:
doubleGyre.flow = set_flow_resolution([2 1]*100,doubleGyre.flow);
doubleGyre.strainline = reset_strainline(doubleGyre.strainline);
[doubleGyre,hAxes] = strain_lcs_script(doubleGyre);
set(findobj(hAxes,'tag','strainline'),'visible','on')

%%
% Furthermore, the strainline integration error tolerance should be
% decreased:
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);
doubleGyre = strain_lcs_script(doubleGyre);
set(findobj(hAxes,'tag','strainline'),'visible','on')

%%
% Now that the strainlines have a good appearance, filtering parameters
% are adjusted to find significant hyperbolic barriers
doubleGyre.strainline = set_strainline_filtering_parameters(struct('distance',1.5,'resolution',[1 1]),doubleGyre.strainline);
doubleGyre = strain_lcs_script(doubleGyre);

%%
% The double gyre is analyzed in
% <http://link.aip.org/link/doi/10.1063/1.3690153 DOI:10.1063/1.3690153>.
% To produce a figure similar to Figure 10, set the timespan to match and
% increase the strainline resolution:
doubleGyre.flow = set_flow_timespan([0 20],doubleGyre.flow);
doubleGyre.strainline = set_strainline_resolution([2 1]*10,doubleGyre.strainline);
doubleGyre = strain_lcs_script(doubleGyre);
