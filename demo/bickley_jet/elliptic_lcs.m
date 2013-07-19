%% Define input parameters
u = 62.66;
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
[cgEigenvalue,cgEigenvector] = eig_cgStrain(bickleyJet.flow,method,customEigMethod,coupledIntegration);

%% Plot finite-time Lyapunov exponent
% FIXME Should have uniform format for cgEigenvalue either m-by-n array or
% two column array
cgEigenvalueArray = reshape(cgEigenvalue,[fliplr(bickleyJet.flow.resolution),2]);
ftle = compute_ftle(cgEigenvalueArray(:,:,2),diff(bickleyJet.flow.timespan));
hAxes = setup_figure(bickleyJet.flow.domain);
plot_ftle(hAxes,bickleyJet.flow,ftle);
drawnow

%% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection{1}.endPosition = [6.5,-1.4;4,-2]*1e6;
poincareSection{2}.endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection{3}.endPosition = [3.25,1.5;1.4,2.6]*1e6;
poincareSection{4}.endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection{5}.endPosition = [1.65e7,1.5e6;1.5e7,2.6e6];
rOrbit = nan(1,2);
for idx = 1:numel(poincareSection)
    poincareSection{idx}.numPoints = 80; %#ok<SAGROW>
    % set maximum integration length to twice the expected circumference
    rOrbit = norm(poincareSection{idx}.endPosition(2,:)-poincareSection{idx}.endPosition(1,:));
    poincareSection{idx}.integrationLength = [0,2*(2*pi*rOrbit)]; %#ok<SAGROW>
end

%% Plot Poincare sections
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),1:numel(poincareSection));
set(hPoincareSection,'color','w')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'markerFaceColor','w')
drawnow

%% Find closed orbits with Poincare section return map
[shearline.etaPos,shearline.etaNeg] = lagrangian_shear(cgEigenvector,cgEigenvalue);
showGraph = true;
odeSolverOptions = odeset('relTol',1e-4);
nBisection = 3;
dThresh = 1e-2;
closedOrbits = poincare_closed_orbit_multi(bickleyJet.flow,shearline,poincareSection,odeSolverOptions,dThresh,showGraph);
% etaPos closed orbits
hClosedOrbitPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}(:,1),closedOrbits{i}{1}(:,2)),1:numel(closedOrbits));
set(hClosedOrbitPos,'color','w')
set(hClosedOrbitPos,'linewidth',2)
% etaNeg closed orbits
hClosedOrbitNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}(:,1),closedOrbits{i}{2}(:,2)),1:numel(closedOrbits));
set(hClosedOrbitNeg,'color','w')
set(hClosedOrbitNeg,'linewidth',2)
