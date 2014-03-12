u = 62.66;

lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];

flow.imposeIncompressibility = true;
flow.periodicBc = [true,false];
perturbationCase = 3;
phiTimespan = [0,25];
phiInitial = [0,0];
phiSol = ode45(@d_phi,phiTimespan,phiInitial);
timeResolution = 1e5;
phi1 = deval(phiSol,linspace(phiTimespan(1),phiTimespan(2),timeResolution),1);
phi1Max = max(phi1);
flow.derivative = @(t,x,~)derivative(t,x,false,u,lengthX,lengthY,epsilon,perturbationCase,phiSol,phi1Max);

flow.domain = [0,lengthX;[-1,1]*2.2599*lengthY];
flow.timespan = [0,4*lengthX/u];
% Make grid Cartesian
resolutionX = 20;
gridSpace = diff(flow.domain(1,:))/(double(resolutionX)-1);
resolutionY = round(diff(flow.domain(2,:))/gridSpace);
flow.resolution = [resolutionX,resolutionY];

flow = animate_flow(flow);
