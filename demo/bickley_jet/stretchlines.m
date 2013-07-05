%% Input parameters
u = 62.66;
earthRadius = 6.371e6;

lengthX = pi*earthRadius;
lengthY = 1.77e6;

bickleyJet.flow.imposeIncompressibility = true;
bickleyJet.flow.periodicBc = [true,false];
perturbationCase = 3;
bickleyJet.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,lengthX,lengthY,perturbationCase),bickleyJet.flow);
bickleyJet.flow = set_flow_resolution([500,200],bickleyJet.flow);
% magicNumber gives the domain an aspect ratio similar to that used in
% doi:10.1016/j.physd.2012.06.012 and ensures grid spacing is equal in the
% x and y directions. In doi:10.1016/j.physd.2012.06.012, 
% magicNumber = 2.2599.
magicNumber = .5*pi*earthRadius/lengthY*double(bickleyJet.flow.resolution(2)-1)/double(bickleyJet.flow.resolution(1)-1);
bickleyJet.flow = set_flow_domain([0,lengthX;[-1,1]*magicNumber*lengthY],bickleyJet.flow);
bickleyJet.flow = set_flow_timespan([0,4*lengthX/u],bickleyJet.flow);

%% Compute Cauchy-Green strain eigenvalues and eigenvectors
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[bickleyJet.flow.cgEigenvalue,bickleyJet.flow.cgEigenvector] = eig_cgStrain(bickleyJet.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue = reshape(bickleyJet.flow.cgEigenvalue,[fliplr(bickleyJet.flow.resolution),2]);
cgEigenvector = reshape(bickleyJet.flow.cgEigenvector,[fliplr(bickleyJet.flow.resolution),4]);

% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(bickleyJet.flow.timespan));
hAxes = setup_figure(bickleyJet.flow.domain);
hImagesc = imagesc(bickleyJet.flow.domain(1,:),bickleyJet.flow.domain(2,:),ftle);
set(hImagesc,'parent',hAxes)
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
drawnow

%% Compute stretchlines
bickleyJet.stretchline.maxLength = 1e7;
bickleyJet.stretchline.resolution = uint64([20,10]);
bickleyJet.stretchline.odeSolverOptions = odeset('RelTol',1e-6);
bickleyJet.stretchline.position = compute_stretchline(bickleyJet.flow,bickleyJet.stretchline);
% Plot stretchlines
hStretchline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),bickleyJet.stretchline.position);
grayColor = [.7,.7,.7];
set(hStretchline,'color',grayColor)
% Calculate relative stretching
segmentIndex = cellfun(@(position)[1,size(position,1)],bickleyJet.stretchline.position,'UniformOutput',false);
gridPosition = initialize_ic_grid(bickleyJet.flow.resolution,bickleyJet.flow.domain);
relativeStretching = relative_stretching(bickleyJet.stretchline.position,segmentIndex,bickleyJet.flow.cgEigenvalue(:,2),bickleyJet.flow.domain,bickleyJet.flow.resolution,false);
relativeStretching = cell2mat(relativeStretching);

[~,sortIndex] = sort(relativeStretching);

% Highlight most-stretching stretchlines
nMost = 50;
hStretchlineMost = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),bickleyJet.stretchline.position(sortIndex(end-nMost:end)));
set(hStretchlineMost,'color','r')
