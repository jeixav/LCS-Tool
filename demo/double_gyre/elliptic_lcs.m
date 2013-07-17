%% Define input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);

doubleGyre.flow = set_flow_domain([0,2;0,1],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([501,251],doubleGyre.flow);

%% Compute Cauchy-Green strain eigenvalues and eigenvectors
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);

%% Plot finite-time Lyapunov exponent
% FIXME Should have uniform format for cgEigenvalue either m-by-n array or
% two column array
cgEigenvalueArray = reshape(cgEigenvalue,[fliplr(doubleGyre.flow.resolution),2]);
ftle = compute_ftle(cgEigenvalueArray(:,:,2),diff(doubleGyre.flow.timespan));
hAxes = setup_figure(doubleGyre.flow.domain);
plot_ftle(hAxes,doubleGyre.flow,ftle);
drawnow

%% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection{1}.endPosition = [.5,.6;.35,.5];
poincareSection{1}.numPoints = 80;
% set maximum integration length to twice the expected circumference
rOrbit = hypot(diff(poincareSection{1}.endPosition(:,1)),diff(poincareSection{1}.endPosition(:,2)));
poincareSection{1}.integrationLength = [0,2*(2*pi*rOrbit)];

%% Plot Poincare sections
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),numel(poincareSection));
set(hPoincareSection,'color','w')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'markerFaceColor','w')
drawnow

%% Find closed orbits with Poincare section return map
[etaPos,etaNeg] = lagrangian_shear(cgEigenvector,cgEigenvalue);
showGraph = true;
odeSolverOptions = odeset('relTol',1e-3);
nBisection = 2;
dThresh = 1e-2;
closedOrbitPos = poincare_closed_orbit(doubleGyre.flow,etaPos,poincareSection{1},odeSolverOptions,nBisection,dThresh,showGraph);
hClosedOrbit = plot(hAxes,closedOrbitPos(:,1),closedOrbitPos(:,2));
set(hClosedOrbit,'color','w')
set(hClosedOrbit,'linewidth',2)

% Repeat for second Poincare section, with etaNeg vector
poincareSection{1}.endPosition = [1.5,.4;1.7,.5];
rOrbit = hypot(diff(poincareSection{1}.endPosition(:,1)),diff(poincareSection{1}.endPosition(:,2)));
poincareSection{1}.integrationLength = [0,2*(2*pi*rOrbit)];
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),numel(poincareSection));
set(hPoincareSection,'color','w')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'markerFaceColor','w')
drawnow
closedOrbitPos = poincare_closed_orbit(doubleGyre.flow,etaNeg,poincareSection{1},odeSolverOptions,nBisection,dThresh,showGraph);
hClosedOrbit = plot(hAxes,closedOrbitPos(:,1),closedOrbitPos(:,2));
set(hClosedOrbit,'color','w')
set(hClosedOrbit,'linewidth',2)
