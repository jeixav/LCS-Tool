function doubleGyre = double_gyre
% Describes a double gyre flow.
%
% References: doi:10.1016/j.physd.2005.10.007, doi:10.5194/npg-7-59-2000,
% doi:10.5194/npg-4-223-1997

t = sym('t');
x = sym('x');
y = sym('y');

p = struct('epsilon',.1,'a',.1,'omega',pi/5);

forcing = p.epsilon*sin(p.omega*t)*x^2 + (1 - 2*p.epsilon...
    *sin(p.omega*t))*x;

doubleGyre.flow.symDerivative(1) = -pi*p.a*sin(pi*forcing)*cos(pi*y);
doubleGyre.flow.symDerivative(2) = pi*p.a*cos(pi*forcing).*sin(pi*y)...
    *(2*p.epsilon*sin(p.omega*t)*x + 1 - 2*p.epsilon*sin(p.omega*t));

doubleGyre.flow = set_flow_domain([0 2; 0 1],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0 20],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([2 1]*10,doubleGyre.flow);

doubleGyre.flow.isCompressible = false;

doubleGyre.strainline = set_strainline_resolution([2 1]*5);
doubleGyre.strainline = set_strainline_max_length(10,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_geodesic_deviation_tol(inf,...
    doubleGyre.strainline);
doubleGyre.strainline = set_strainline_length_tol(0,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_filtering_method('superminimize',...
    doubleGyre.strainline);
doubleGyre.strainline = set_strainline_filtering_parameters(...
    struct('distance',0,'resolution',[1 1]),doubleGyre.strainline);

doubleGyre.shearline = set_shearline_resolution(uint64([2 1]*5));
doubleGyre.shearline = set_shearline_max_length(10,doubleGyre.shearline);
doubleGyre.shearline = set_shearline_average_geodesic_deviation_tol(...
    [inf inf],doubleGyre.shearline);
