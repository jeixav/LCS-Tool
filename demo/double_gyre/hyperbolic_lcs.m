epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)double_gyre_derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);

doubleGyre.flow = set_flow_domain([0,2;0,1],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([2,1]*100,doubleGyre.flow);

doubleGyre.strainline = set_strainline_resolution([2,1]*5);
doubleGyre.strainline = set_strainline_max_length(5,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);
doubleGyre.strainline = set_strainline_geodesic_deviation_tol(inf,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_length_tol(0,doubleGyre.strainline);
doubleGyre.strainline = set_strainline_filtering_method('superminimize',doubleGyre.strainline);
doubleGyre.strainline = set_strainline_filtering_parameters(struct('distance',1.5,'resolution',[1,1]),doubleGyre.strainline);

doubleGyre = strain_lcs_script(doubleGyre);
