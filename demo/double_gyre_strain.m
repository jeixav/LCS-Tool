function [doubleGyre,hStrainAxes] = double_gyre_strain

doubleGyre = double_gyre;

doubleGyre.flow = set_flow_resolution(500,doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
doubleGyre.flow = set_flow_ode_solver_options(odeset('RelTol',1e-8,'AbsTol',1e-8),doubleGyre.flow);
doubleGyre.flow.cgStrainMethod.name = 'finiteDifference';
doubleGyre.flow.coupledIntegration = 1e4;

doubleGyre.strainline = set_strainline_filtering_parameters(struct('distance',3,'resolution',[2,2]),doubleGyre.strainline);

showPlot.strainline = true;

verbose.progress = false;

[doubleGyre,hStrainAxes] = strain_lcs_script(doubleGyre,showPlot,verbose);
