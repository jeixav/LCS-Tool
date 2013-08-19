u = 62.66;

lengthX = pi*earthRadius;
lengthY = 1.77e6;

flow.imposeIncompressibility = true;
flow.periodicBc = [true,false];
perturbationCase = 3;
flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,lengthX,lengthY,perturbationCase),flow);

flow = set_flow_domain([0,lengthX;[-1,1]*2.2599*lengthY],flow);
flow = set_flow_timespan([0,4*lengthX/u],flow);
flow = set_flow_resolution(20,flow);

flow = animate_flow(flow);
