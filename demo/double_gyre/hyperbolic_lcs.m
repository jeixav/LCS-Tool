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

doubleGyre.strainline = set_strainline_max_length(20);
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);

gridSpace = diff(doubleGyre.flow.domain(1,:))/(double(doubleGyre.flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;

%% Compute Cauchy-Green strain eigenvalues and eigenvectors
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue = reshape(cgEigenvalue,[fliplr(doubleGyre.flow.resolution),2]);
cgEigenvector = reshape(cgEigenvector,[fliplr(doubleGyre.flow.resolution),4]);

%% Repelling LCS analysis
% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(doubleGyre.flow.timespan));
hAxes = setup_figure(doubleGyre.flow.domain);
hImagesc = imagesc(doubleGyre.flow.domain(1,:),doubleGyre.flow.domain(2,:),ftle);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow
% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyre.strainline.maxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),doubleGyre.flow.domain);
% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','w')
set(hStrainline,'lineWidth',2)
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:numel(strainlinePosition));
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','k')
set(hStrainlineInitialPosition,'MarkerFaceColor','k')

%% Attracting LCS analysis
% Plot backward time finite-time Lyapunov exponent
doubleGyreBackward = doubleGyre;
doubleGyreBackward.flow = set_flow_timespan([0,-timespan],doubleGyre.flow);
cgEigenvalueBackward = eig_cgStrain(doubleGyreBackward.flow,method,customEigMethod,coupledIntegration);
cgEigenvalueBackward2 = reshape(cgEigenvalueBackward(:,2),fliplr(doubleGyreBackward.flow.resolution));
ftleBackward = compute_ftle(cgEigenvalueBackward2,diff(doubleGyreBackward.flow.timespan));
hAxes = setup_figure(doubleGyreBackward.flow.domain);
hImagesc = imagesc(doubleGyreBackward.flow.domain(1,:),doubleGyreBackward.flow.domain(2,:),ftleBackward);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow
% Compute stretchlines
nMaxStretchlines = uint8(40);
stretchlineMaxLength = doubleGyre.strainline.maxLength;
[stretchlinePosition,stretchlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,:,1),cgEigenvector(:,:,3:4),doubleGyre.flow.domain,nMaxStretchlines);
% Plot stretchlines
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlinePosition);
set(hStretchline,'color','w')
set(hStretchline,'lineWidth',2)
hStretchlineInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineInitialPosition(1,idx),stretchlineInitialPosition(2,idx)),1:numel(stretchlinePosition));
set(hStretchlineInitialPosition,'marker','o')
set(hStretchlineInitialPosition,'MarkerEdgeColor','k')
set(hStretchlineInitialPosition,'MarkerFaceColor','k')
