epsilon = .1;
amplitude = .1;
omega = pi/5;
flow.imposeIncompressibility = true;
flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),flow);

flow = set_flow_domain([0,2;0,1],flow);
flow = set_flow_timespan([0,20],flow);
flow = set_flow_resolution([2,1]*10,flow);

flow = animate_flow(flow);
