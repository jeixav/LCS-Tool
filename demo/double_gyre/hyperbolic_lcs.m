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

%% Repelling LCS analysis
% Compute Cauchy-Green strain eigenvalues and eigenvectors
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue = reshape(cgEigenvalue,[fliplr(doubleGyre.flow.resolution),2]);
cgEigenvector = reshape(cgEigenvector,[fliplr(doubleGyre.flow.resolution),4]);

% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(doubleGyre.flow.timespan));
hAxes = setup_figure(doubleGyre.flow.domain);
plot_ftle(hAxes,doubleGyre.flow,ftle);
hImagesc = imagesc(doubleGyre.flow.domain(1,:),doubleGyre.flow.domain(2,:),ftle);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyre.strainline.maxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),doubleGyre.flow.domain,doubleGyre.flow.periodicBc);

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','r')
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:numel(strainlinePosition));
set(hStrainlineInitialPosition,'MarkerSize',2)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor','r')

%% Attracting LCS analysis

% Compute Cauchy-Green strain eigenvalues and eigenvectors
doubleGyreBackward = doubleGyre;
doubleGyreBackward.flow = set_flow_timespan([timespan,0],doubleGyre.flow);
[cgEigenvalueBackward,cgEigenvectorBackward] = eig_cgStrain(doubleGyreBackward.flow,method,customEigMethod,coupledIntegration);
cgEigenvalueBackward = reshape(cgEigenvalueBackward,[fliplr(doubleGyreBackward.flow.resolution),2]);
cgEigenvectorBackward = reshape(cgEigenvectorBackward,[fliplr(doubleGyreBackward.flow.resolution),4]);

% Plot backward time finite-time Lyapunov exponent
ftleBackward = compute_ftle(cgEigenvalueBackward(:,:,2),diff(doubleGyreBackward.flow.timespan));
hAxes = setup_figure(doubleGyreBackward.flow.domain);
plot_ftle(hAxes,doubleGyreBackward.flow,ftleBackward);
drawnow

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyreBackward.strainline.maxLength,cgEigenvalueBackward(:,:,2),cgEigenvectorBackward(:,:,1:2),doubleGyre.flow.domain,doubleGyre.flow.periodicBc);

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','b')
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:numel(strainlinePosition));
set(hStrainlineInitialPosition,'MarkerSize',2)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor','b')
