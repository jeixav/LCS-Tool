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
hImagesc = imagesc(doubleGyre.flow.domain(1,:),doubleGyre.flow.domain(2,:),ftle);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow

%% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection{1}.endPosition = [.5,.6;.25,.5];
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
odeSolverOptions = [];
[closedOrbitPos,orbitsPos] = poincare_closed_orbit(doubleGyre.flow,etaPos,poincareSection{1},[],showGraph);
hClosedOrbit = plot(hAxes,closedOrbitPos(:,1),closedOrbitPos(:,2));
set(hClosedOrbit,'color','w')
set(hClosedOrbit,'linewidth',2)
