epsilon = .1;
amplitude = .1;
omega = pi/5;
flow.imposeIncompressibility = true;
flow.derivative = @(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega);

flow.domain = [0,2;0,1];
flow.timespan = [0,10];
flow.resolution = [2,1]*20;

animate_flow(flow);
