% duffing Defines a Duffing oscillator
%
% EXAMPLE
% matlabpool('open')
% pctRunOnAll javaaddpath('ParforProgress2')
% Duffing = duffing;
% showPlot.strainline = true;
% Duffing = strain_lcs_script(Duffing,showPlot);

function Duffing = duffing

t = sym('t');
x = sym('x');
y = sym('y');

flow.symDerivative(1) = y;
flow.symDerivative(2) = x - x^3 - .15*y + .3*cos(t);

flow = set_flow_domain([-2 2;-1 1.5],flow);
flow = set_flow_timespan([8*pi 0],flow);
flow = set_flow_resolution(200,flow);

strainline = set_strainline_resolution([2 1]*4);
strainline = set_strainline_max_length(20,strainline);

strainline = set_strainline_geodesic_deviation_tol(inf,strainline);
strainline = set_strainline_length_tol(0,strainline);
strainline = set_strainline_ode_solver_options(odeset('relTol',1e-4),...
    strainline);
strainline = set_strainline_filtering_method('superminimize',strainline);
filteringParameters.distance = 0;
filteringParameters.resolution = [2 2];
strainline = set_strainline_filtering_parameters(filteringParameters,...
    strainline);

Duffing = struct('flow',flow,'strainline',strainline);
