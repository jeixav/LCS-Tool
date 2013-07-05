%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);
timespan = 20;

doubleGyre.flow = set_flow_domain([-.1,2.1;-.05,1.05],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,timespan],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([1101,551],doubleGyre.flow);

doubleGyre.strainline = set_strainline_max_length(20);
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);

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
doubleGyre.stretchline.position = compute_stretchline(doubleGyre.flow,doubleGyre.stretchline);
% Plot stretchlines
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),doubleGyre.stretchline.position);
grayColor = [.7,.7,.7];
set(hStretchline,'color',grayColor)
% Calculate relative stretching
segmentIndex = cellfun(@(position)[1,size(position,1)],doubleGyre.stretchline.position,'UniformOutput',false);
gridPosition = initialize_ic_grid(doubleGyre.flow.resolution,doubleGyre.flow.domain);
relativeStretching = relative_stretching(doubleGyre.stretchline.position,segmentIndex,doubleGyre.flow.cgEigenvalue(:,2),doubleGyre.flow.domain,doubleGyre.flow.resolution,false);
relativeStretching = cell2mat(relativeStretching);

[~,sortIndex] = sort(relativeStretching);

% Highlight most-stretching stretchlines
nMost = 20;
hStretchlineMost = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),doubleGyre.stretchline.position(sortIndex(end-nMost:end)));
set(hStretchlineMost,'color','r')
