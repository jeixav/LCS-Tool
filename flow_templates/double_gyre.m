% double_gyre Defines double gyre flow for LCS analysis
%
% REFERENCES
% doi:10.1016/j.physd.2005.10.007
% doi:10.5194/npg-7-59-2000,
% doi:10.5194/npg-4-223-1997
%
% EXAMPLE
% matlabpool('open')
% pctRunOnAll javaaddpath('ParforProgress2')
% addpath('flow_templates')
% doubleGyre = double_gyre;
% doubleGyre = strain_lcs_script(doubleGyre);
% set(findobj(gca,'tag','strainline'),'visible','on')

function doubleGyre = double_gyre

epsilon = .07;
amplitude = .1;
omega = pi/5;

doubleGyre.flow.derivative = @(t,x,useEoV)double_gyre_derivative(t,x,useEoV,epsilon,amplitude,omega);
doubleGyre.flow.cgStrainMethod.name = 'equationOfVariation';

doubleGyre.flow = set_flow_domain([0,2;0,1],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([2,1]*100,doubleGyre.flow);

doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow.coupledIntegration = 1e5;

doubleGyre.strainline = set_strainline_resolution([2,1]*10);
doubleGyre.strainline = set_strainline_max_length(10,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_geodesic_deviation_tol(inf,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_length_tol(0,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_filtering_method('superminimize',doubleGyre.strainline);
doubleGyre.strainline = set_strainline_filtering_parameters(struct('distance',1.5,'resolution',[1,1]),doubleGyre.strainline);
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);

doubleGyre.shearline = set_shearline_resolution([2,1]*10);
doubleGyre.shearline = set_shearline_max_length(10,doubleGyre.shearline);
doubleGyre.shearline = set_shearline_average_geodesic_deviation_tol([inf inf],doubleGyre.shearline);
doubleGyre.shearline = set_shearline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.shearline);

end