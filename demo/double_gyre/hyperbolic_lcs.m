%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);

doubleGyre.flow = set_flow_domain([-.1,2.1;-.05,1.05],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([551,276],doubleGyre.flow);

doubleGyre.strainline = set_strainline_max_length(10);
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);

gridSpace = diff(doubleGyre.flow.domain(1,:))/(double(doubleGyre.flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;

%% Compute λ₂ and ξ₁
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(doubleGyre.flow.resolution));
cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(doubleGyre.flow.resolution),2]);

%% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue2,diff(doubleGyre.flow.timespan));
hAxes = setup_figure(doubleGyre.flow.domain);
hImagesc = imagesc(doubleGyre.flow.domain(1,:),doubleGyre.flow.domain(2,:),ftle);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow

%% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda2_max(localMaxDistance,doubleGyre.strainline.maxLength,cgEigenvalue2,cgEigenvector1,doubleGyre.flow.domain);

%% Plot strainlines
hStrainline = cellfun(@(position)plot(position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','w')
set(hStrainline,'lineWidth',2)
hStrainlineInitialPosition = arrayfun(@(idx)plot(strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:numel(strainlinePosition));
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','k')
set(hStrainlineInitialPosition,'MarkerFaceColor','k')
