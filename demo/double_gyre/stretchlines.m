%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);
timespan = 20;
doubleGyre.flow = set_flow_domain([-.1,2.1;-.05,1.05],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([551,276],doubleGyre.flow);

doubleGyre.stretchline.maxLength = 2;
doubleGyre.stretchline.resolution = uint64([20,10]);
doubleGyre.stretchline.odeSolverOptions = odeset('RelTol',1e-6);

% Cauchy-Green strain numerical parameters
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;

% Plot parameters
forwardLcsColor = 'b';
backwardLcsColor = 'r';

%% Forward-time stretchlines
doubleGyre.flow = set_flow_timespan([0,timespan],doubleGyre.flow);
[doubleGyre.flow.cgEigenvalue,doubleGyre.flow.cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
doubleGyre.stretchline.position = compute_stretchline(doubleGyre.flow,doubleGyre.stretchline);
hAxes = setup_figure(doubleGyre.flow.domain);
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),doubleGyre.stretchline.position);
set(hStretchline,'color',forwardLcsColor)
title(hAxes,'Forward-time Stretchlines')
drawnow

%% Backward-time stretchlines
doubleGyre.flow = set_flow_timespan([timespan,0],doubleGyre.flow);
[doubleGyre.flow.cgEigenvalue,doubleGyre.flow.cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
doubleGyre.stretchline.position = compute_stretchline(doubleGyre.flow,doubleGyre.stretchline);
hAxes = setup_figure(doubleGyre.flow.domain);
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),doubleGyre.stretchline.position);
set(hStretchline,'color',backwardLcsColor)
title(hAxes,'Backward-time Stretchlines')
