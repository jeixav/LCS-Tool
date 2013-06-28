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

%% Compute stretchlines
nMaxStretchlines = uint8(40);
stretchlineMaxLength = doubleGyre.strainline.maxLength;
[stretchlinePosition,stretchlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,:,1),cgEigenvector(:,:,3:4),doubleGyre.flow.domain,nMaxStretchlines);
% Plot stretchlines
hAxes = setup_figure(doubleGyre.flow.domain);
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlinePosition);
set(hStretchline,'color','k')
hStretchlineInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineInitialPosition(1,idx),stretchlineInitialPosition(2,idx)),1:numel(stretchlinePosition));
set(hStretchlineInitialPosition,'marker','o')
set(hStretchlineInitialPosition,'MarkerEdgeColor','k')
set(hStretchlineInitialPosition,'MarkerFaceColor','k')
% Advect stretchlines under forward-time flow map
flowSolution = cellfun(@(position)integrate_flow(doubleGyre.flow,position),stretchlinePosition,'uniformOutput',false);
advectedPosition = cellfun(@(solution)deval(solution,doubleGyre.flow.timespan(end)),flowSolution,'uniformOutput',false);
advectedPosition = cellfun(@(position)reshape(position,2,numel(position)/2),advectedPosition,'uniformOutput',false);
% Plot advected stretchlines
hStretchlineAdv = cellfun(@(position)plot(hAxes,position(1,:),position(2,:)),advectedPosition);
set(hStretchlineAdv,'color','k')
set(hStretchlineAdv,'linestyle','none')
set(hStretchlineAdv,'marker','.')
