%% Input parameters
u = 62.66;

lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];

bickleyJet.flow.imposeIncompressibility = true;
bickleyJet.flow.periodicBc = [true,false];
perturbationCase = 3;
bickleyJet.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,u,lengthX,lengthY,epsilon,perturbationCase),bickleyJet.flow);
bickleyJet.flow = set_flow_resolution([500,200],bickleyJet.flow);
% magicNumber gives the domain an aspect ratio similar to that used in
% doi:10.1016/j.physd.2012.06.012 and ensures grid spacing is equal in the
% x and y directions. In doi:10.1016/j.physd.2012.06.012, 
% magicNumber = 2.2599.
magicNumber = .5*pi*earthRadius/lengthY*double(bickleyJet.flow.resolution(2)-1)/double(bickleyJet.flow.resolution(1)-1);
bickleyJet.flow = set_flow_domain([0,lengthX;[-1,1]*magicNumber*lengthY],bickleyJet.flow);
timespan = 4*lengthX/u;

bickleyJet.stretchline.maxLength = 1e7;
bickleyJet.stretchline.resolution = uint64([20,10]);
bickleyJet.stretchline.odeSolverOptions = odeset('RelTol',1e-6);

% Cauchy-Green strain numerical parameters
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;

% Plot parameters
forwardLcsColor = 'b';
backwardLcsColor = 'r';

%% Forward-time stretchlines
bickleyJet.flow = set_flow_timespan([0,4*lengthX/u],bickleyJet.flow);
[bickleyJet.flow.cgEigenvalue,bickleyJet.flow.cgEigenvector] = eig_cgStrain(bickleyJet.flow,method,customEigMethod,coupledIntegration);
bickleyJet.stretchline.position = compute_stretchline(bickleyJet.flow,bickleyJet.stretchline);
hAxes = setup_figure(bickleyJet.flow.domain);
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),bickleyJet.stretchline.position);
set(hStretchline,'color',forwardLcsColor)
title(hAxes,'Forward-time Stretchlines')
drawnow

%% Backward-time stretchlines
bickleyJet.flow = set_flow_timespan([4*lengthX/u,0],bickleyJet.flow);
[bickleyJet.flow.cgEigenvalue,bickleyJet.flow.cgEigenvector] = eig_cgStrain(bickleyJet.flow,method,customEigMethod,coupledIntegration);
bickleyJet.stretchline.position = compute_stretchline(bickleyJet.flow,bickleyJet.stretchline);
hAxes = setup_figure(bickleyJet.flow.domain);
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),bickleyJet.stretchline.position);
set(hStretchline,'color',backwardLcsColor)
title(hAxes,'Backward-time Stretchlines')
