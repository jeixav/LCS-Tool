%% Lagrangian Coherent Structures Toolbox User Guide
%

%% Introduction
% Use of the LCS toolbox is demonstrated by analyzing a double gyre flow:
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
% The flow vector field is defined as a function. Here is the definition of
% a double gyre:
t = sym('t');
x = sym('x');
y = sym('y');

p = struct('epsilon',.1,'a',.1,'omega',pi/5);

forcing = p.epsilon*sin(p.omega*t)*x^2 + (1 - 2*p.epsilon...
    *sin(p.omega*t))*x;

flow.symDerivative(1) = -pi*p.a*sin(pi*forcing)*cos(pi*y);
flow.symDerivative(2) = pi*p.a*cos(pi*forcing).*sin(pi*y)...
    *(2*p.epsilon*sin(p.omega*t)*x + 1 - 2*p.epsilon*sin(p.omega*t));

flow.derivative = sym2fun(flow.symDerivative);

%%
% The flow domain, timespan and resolution must also be defined:
flow.domain = [0 2; 0 1];
flow.timespan = [0 20];
flow.resolution = uint64([2 1]*10);

%% Flow animation
% To verify that the flow has been correctly defined, it can be animated:
animate_flow(flow);

%%
% Parameters can be changed and the animation re-run
% ...

%% Hyperbolic barriers
%

strainline.resolution = uint64([2 1]*5);
strainline.finalTime = 10;
strainline.geodesicDeviationTol = inf;
strainline.lengthTol = 0;
strainline.filteringMethod = 'hausdorff';
strainline.filteringDistanceTol = 0;
strainline.odeSolverOptions = odeset('relTol',1e-4);
doubleGyre = struct('flow',flow,'strainline',strainline);
doubleGyre = strain_lcs_script(doubleGyre);

%% Parabolic barriers
% 

%% Elliptic barriers
%
