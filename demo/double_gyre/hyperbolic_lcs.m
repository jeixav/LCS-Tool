%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);

doubleGyre.flow = set_flow_domain([-.1,2.1;-.05,1.05],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
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
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyre.strainline.maxLength,-cgEigenvalue(:,:,1),cgEigenvector(:,:,1:2),doubleGyre.flow.domain);
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
doubleGyreBackward.flow = set_flow_timespan([20,0],doubleGyre.flow);
cgEigenvalueBackward = eig_cgStrain(doubleGyreBackward.flow,method,customEigMethod,coupledIntegration);
cgEigenvalueBackward2 = reshape(cgEigenvalueBackward(:,2),fliplr(doubleGyreBackward.flow.resolution));
ftleBackward = compute_ftle(cgEigenvalueBackward2,diff(doubleGyreBackward.flow.timespan));
hAxes = setup_figure(doubleGyreBackward.flow.domain);
hImagesc = imagesc(doubleGyreBackward.flow.domain(1,:),doubleGyreBackward.flow.domain(2,:),ftleBackward);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow
% Low-resolution quiver plot of ξ₂
doubleGyreLowRes = doubleGyre;
doubleGyreLowRes.flow = set_flow_resolution([101,51],doubleGyre.flow);
[~,cgEigenvectorLowRes] = eig_cgStrain(doubleGyreLowRes.flow,method,customEigMethod,coupledIntegration);
cgEigenvectorLowRes2 = reshape(cgEigenvectorLowRes(:,3:4),[fliplr(doubleGyreLowRes.flow.resolution),2]);
hAxesXi2 = setup_figure(doubleGyreLowRes.flow.domain);
gridPosition{1} = linspace(doubleGyreLowRes.flow.domain(1,1),doubleGyreLowRes.flow.domain(1,2),doubleGyreLowRes.flow.resolution(1));
gridPosition{2} = linspace(doubleGyreLowRes.flow.domain(2,1),doubleGyreLowRes.flow.domain(2,2),doubleGyreLowRes.flow.resolution(2));
hQuiver = quiver(hAxesXi2,gridPosition{1},gridPosition{2},cgEigenvectorLowRes2(:,:,1),cgEigenvectorLowRes2(:,:,2));
set(hQuiver,'AutoScaleFactor',.5)
% Compute strainlines
nMaxStretchlines = uint8(10);
[stretchlinePosition,stretchlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyre.strainline.maxLength,-cgEigenvalue(:,:,1),cgEigenvector(:,:,3:4),doubleGyre.flow.domain,nMaxStretchlines);
% Plot strainlines
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlinePosition);
set(hStretchline,'color','w')
set(hStretchline,'lineWidth',2)
hStretchlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:numel(strainlinePosition));
set(hStretchlineInitialPosition,'marker','o')
set(hStretchlineInitialPosition,'MarkerEdgeColor','k')
set(hStretchlineInitialPosition,'MarkerFaceColor','k')
