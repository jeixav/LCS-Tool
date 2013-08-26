%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);
timespan = 20;

doubleGyre.flow = set_flow_domain([-.1,2.1;-.05,1.05],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,timespan],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([551,276],doubleGyre.flow);
doubleGyre.flow.periodicBc = [false,false];

doubleGyre.strainline = set_strainline_max_length(20);
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);

gridSpace = diff(doubleGyre.flow.domain(1,:))/(double(doubleGyre.flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;

forwardLcsColor = 'r';
backwardLcsColor = 'b';
shearLcsColor = [0,.6,0];

%% Forward-time LCS analysis
% Compute Cauchy-Green strain eigenvalues and eigenvectors
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[doubleGyre.flow.cgEigenvalue,doubleGyre.flow.cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
% FIXME Should use m-by-n or (m*n)-by-2 array forms throughout LCS Tool
cgEigenvalue = reshape(doubleGyre.flow.cgEigenvalue,[fliplr(doubleGyre.flow.resolution),2]);
cgEigenvector = reshape(doubleGyre.flow.cgEigenvector,[fliplr(doubleGyre.flow.resolution),4]);

% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(doubleGyre.flow.timespan));
hAxes = setup_figure(doubleGyre.flow.domain);
title(hAxes,'Forward-time LCS')
plot_ftle(hAxes,doubleGyre.flow,ftle);
colormap(hAxes,flipud(gray))
drawnow

% Compute closed shearlines
% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection{1}.endPosition = [.5,.6;.35,.5];
poincareSection{1}.numPoints = 80;
% set maximum integration length to twice the expected circumference
rOrbit = hypot(diff(poincareSection{1}.endPosition(:,1)),diff(poincareSection{1}.endPosition(:,2)));
poincareSection{1}.integrationLength = [0,2*(2*pi*rOrbit)];

% Plot Poincare section
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),numel(poincareSection));
set(hPoincareSection,'color',shearLcsColor)
set(hPoincareSection,'linestyle','--')
drawnow

% Find closed orbits with Poincare section return map
[etaPos,etaNeg] = lagrangian_shear(doubleGyre.flow.cgEigenvector,doubleGyre.flow.cgEigenvalue);
showGraph = true;
odeSolverOptions = odeset('relTol',1e-3);
nBisection = 2;
dThresh = 1e-2;
closedOrbitPos = poincare_closed_orbit(doubleGyre.flow,etaPos,poincareSection{1},odeSolverOptions,nBisection,dThresh,showGraph);

% Plot closed orbits
hClosedOrbits = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbitPos);
set(hClosedOrbits,'color',shearLcsColor)
hOuterClosedOrbit = plot(hAxes,closedOrbitPos{end}(:,1),closedOrbitPos{end}(:,2));
set(hOuterClosedOrbit,'color',shearLcsColor)
set(hOuterClosedOrbit,'LineWidth',2)

% Repeat for second Poincare section, with etaNeg vector
poincareSection{1}.endPosition = [1.5,.4;1.7,.5];
rOrbit = hypot(diff(poincareSection{1}.endPosition(:,1)),diff(poincareSection{1}.endPosition(:,2)));
poincareSection{1}.integrationLength = [0,2*(2*pi*rOrbit)];
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),numel(poincareSection));
set(hPoincareSection,'color',shearLcsColor)
set(hPoincareSection,'linestyle','--')
drawnow

closedOrbitNeg = poincare_closed_orbit(doubleGyre.flow,etaNeg,poincareSection{1},odeSolverOptions,nBisection,dThresh,showGraph);

hClosedOrbits = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbitNeg);
set(hClosedOrbits,'color',shearLcsColor)
hOuterClosedOrbit = plot(hAxes,closedOrbitNeg{end}(:,1),closedOrbitNeg{end}(:,2));
set(hOuterClosedOrbit,'color',shearLcsColor)
set(hOuterClosedOrbit,'LineWidth',2)

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyre.strainline.maxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),doubleGyre.flow.domain,doubleGyre.flow.periodicBc);
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbitPos{end});
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbitNeg{end});

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',forwardLcsColor)
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'MarkerSize',2)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor',forwardLcsColor)
drawnow

%% Backward-time LCS analysis
% Compute Cauchy-Green strain eigenvalues and eigenvectors
doubleGyreBackward = doubleGyre;
doubleGyreBackward.flow = set_flow_timespan([timespan,0],doubleGyre.flow);
[doubleGyreBackward.flow.cgEigenvalue,doubleGyreBackward.flow.cgEigenvector] = eig_cgStrain(doubleGyreBackward.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue = reshape(doubleGyreBackward.flow.cgEigenvalue,[fliplr(doubleGyreBackward.flow.resolution),2]);
cgEigenvector = reshape(doubleGyreBackward.flow.cgEigenvector,[fliplr(doubleGyreBackward.flow.resolution),4]);

% Plot backward time finite-time Lyapunov exponent
ftleBackward = compute_ftle(cgEigenvalue(:,:,2),diff(doubleGyreBackward.flow.timespan));
hAxes = setup_figure(doubleGyreBackward.flow.domain);
title(hAxes,'Backward-time LCS')
plot_ftle(hAxes,doubleGyreBackward.flow,ftleBackward);
colormap(hAxes,gray)
drawnow

% Compute closed shearlines
% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection{1}.endPosition = [.5,.6;.3,.5];
poincareSection{1}.numPoints = 80;
% set maximum integration length to twice the expected circumference
rOrbit = hypot(diff(poincareSection{1}.endPosition(:,1)),diff(poincareSection{1}.endPosition(:,2)));
poincareSection{1}.integrationLength = [0,2*(2*pi*rOrbit)];

% Plot Poincare section
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),numel(poincareSection));
set(hPoincareSection,'color',shearLcsColor)
set(hPoincareSection,'linestyle','--')
drawnow

% Find closed orbits with Poincare section return map
[etaPos,etaNeg] = lagrangian_shear(doubleGyreBackward.flow.cgEigenvector,doubleGyreBackward.flow.cgEigenvalue);
showGraph = true;
odeSolverOptions = odeset('relTol',1e-3);
nBisection = 2;
dThresh = 1e-2;
closedOrbitPos = poincare_closed_orbit(doubleGyreBackward.flow,etaNeg,poincareSection{1},odeSolverOptions,nBisection,dThresh,showGraph);

hClosedOrbits = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbitPos);
set(hClosedOrbits,'color',shearLcsColor)
hOuterClosedOrbit = plot(hAxes,closedOrbitPos{end}(:,1),closedOrbitPos{end}(:,2));
set(hOuterClosedOrbit,'color',shearLcsColor)
set(hOuterClosedOrbit,'LineWidth',2)

% Repeat for second Poincare section, with etaPos vector
poincareSection{1}.endPosition = [1.5,.4;1.7,.5];
rOrbit = hypot(diff(poincareSection{1}.endPosition(:,1)),diff(poincareSection{1}.endPosition(:,2)));
poincareSection{1}.integrationLength = [0,2*(2*pi*rOrbit)];
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),numel(poincareSection));
set(hPoincareSection,'color',shearLcsColor)
set(hPoincareSection,'linestyle','--')
drawnow

closedOrbitNeg = poincare_closed_orbit(doubleGyre.flow,etaPos,poincareSection{1},odeSolverOptions,nBisection,dThresh,showGraph);

hClosedOrbits = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbitNeg);
set(hClosedOrbits,'color',shearLcsColor)
hOuterClosedOrbit = plot(hAxes,closedOrbitNeg{end}(:,1),closedOrbitNeg{end}(:,2));
set(hOuterClosedOrbit,'color',shearLcsColor)
set(hOuterClosedOrbit,'LineWidth',2)

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyreBackward.strainline.maxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),doubleGyreBackward.flow.domain,doubleGyreBackward.flow.periodicBc);
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbitPos{end});
strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbitNeg{end});

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',backwardLcsColor)
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'MarkerSize',2)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor',backwardLcsColor)
