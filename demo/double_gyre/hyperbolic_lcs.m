%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
doubleGyre.flow.imposeIncompressibility = true;
doubleGyre.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega),doubleGyre.flow);

doubleGyre.flow = set_flow_domain([0,2;0,1],doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
doubleGyre.flow = set_flow_resolution([501,251],doubleGyre.flow);

doubleGyre.strainline = set_strainline_max_length(5);
doubleGyre.strainline = set_strainline_ode_solver_options(odeset('relTol',1e-6),doubleGyre.strainline);

localMaxDistance = 4/double(doubleGyre.flow.resolution(1)-1);

%% Compute λ₂ and ξ₁
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(doubleGyre.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(doubleGyre.flow.resolution));
cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(doubleGyre.flow.resolution),2]);

%% Plot λ₂ field
hFigure = figure;
hAxes = axes;
set(hAxes,'parent',hFigure)
set(hAxes,'NextPlot','add')
set(hAxes,'box','on')
set(hAxes,'DataAspectRatio',[1,1,1])
hImagesc = imagesc(doubleGyre.flow.domain(1,:),doubleGyre.flow.domain(2,:),log(cgEigenvalue2));
set(hImagesc,'parent',hAxes)
axis(hAxes,'tight')
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','log(\lambda_2)')
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
