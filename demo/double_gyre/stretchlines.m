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
[doubleGyre.flow.cgEigenvalue,doubleGyre.flow.cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue = reshape(doubleGyre.flow.cgEigenvalue,[fliplr(doubleGyre.flow.resolution),2]);
cgEigenvector = reshape(doubleGyre.flow.cgEigenvector,[fliplr(doubleGyre.flow.resolution),4]);

% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(doubleGyre.flow.timespan));
hAxes = setup_figure(doubleGyre.flow.domain);
hImagesc = imagesc(doubleGyre.flow.domain(1,:),doubleGyre.flow.domain(2,:),ftle);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow

%% Compute stretchlines
doubleGyre.stretchline.maxLength = 2;
doubleGyre.stretchline.resolution = uint64([20,10]);
doubleGyre.stretchline.odeSolverOptions = odeset('RelTol',1e-6);
[doubleGyre.stretchline.position,doubleGyre.stretchline.initialPosition] = seed_curves_from_lambda_max(localMaxDistance,doubleGyre.stretchline.maxLength,-cgEigenvalue(:,:,2),cgEigenvector(:,:,3:4),doubleGyre.flow.domain);
% Plot stretchlines
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),doubleGyre.stretchline.position);
set(hStretchline,'color','r')
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,doubleGyre.stretchline.initialPosition(1,idx),doubleGyre.stretchline.initialPosition(2,idx)),1:numel(doubleGyre.stretchline.position));
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','k')
set(hStrainlineInitialPosition,'MarkerFaceColor','k')
