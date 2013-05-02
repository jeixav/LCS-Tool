function [doubleGyre,hStrainAxes] = double_gyre_strain

doubleGyre = double_gyre;

%% xi_1 quiver plot at low resolution
doubleGyre.flow = set_flow_resolution(40,doubleGyre.flow);
doubleGyre.flow = set_flow_ode_solver_options(odeset('RelTol',1e-3,'AbsTol',1e-6),doubleGyre.flow);
doubleGyre.flow.cgStrainMethod.name = 'finiteDifference';
doubleGyre.flow.coupledIntegration = 1e4;

doubleGyre.strainline = set_strainline_max_length(.1,doubleGyre.strainline);

showPlot.quiver = true;
showPlot.strainlineFiltered = false;
verbose.progress = false;

strain_lcs_script(doubleGyre,showPlot,verbose);

%% Strainline plot at high resolution
doubleGyre.flow = set_flow_resolution(500,doubleGyre.flow);
doubleGyre.flow = set_flow_ode_solver_options(odeset('RelTol',1e-8,'AbsTol',1e-8),doubleGyre.flow);
doubleGyre.flow.cgStrainMethod.name = 'finiteDifference';
doubleGyre.flow.coupledIntegration = 1e4;

doubleGyre.strainline = set_strainline_max_length(10,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_filtering_parameters(struct('distance',3,'resolution',[2,2]),doubleGyre.strainline);

showPlot.strainline = true;

verbose.progress = false;

[doubleGyre,hStrainAxes] = strain_lcs_script(doubleGyre,showPlot,verbose);
