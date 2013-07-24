%% Input parameters
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

bickleyJet.strainline = set_strainline_max_length(1e8);
gridSpace = diff(bickleyJet.flow.domain(1,:))/(double(bickleyJet.flow.resolution(1))-1);
localMaxDistance = 8*gridSpace;

%% Repelling LCS analysis
% Compute λ₂ and ξ₁
method.name = 'finiteDifference';
customEigMethod = false;
coupledIntegration = true;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(bickleyJet.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(bickleyJet.flow.resolution));
cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(bickleyJet.flow.resolution),2]);

% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue2,diff(bickleyJet.flow.timespan));
hAxes = setup_figure(bickleyJet.flow.domain);
title(hAxes,'Repelling LCS')
plot_ftle(hAxes,bickleyJet.flow,ftle);
drawnow

% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection{1}.endPosition = [6.5,-1.4;4,-2]*1e6;
poincareSection{2}.endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection{3}.endPosition = [3.25,1.5;1.4,2.6]*1e6;
poincareSection{4}.endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection{5}.endPosition = [1.65e7,1.5e6;1.5e7,2.6e6];
for idx = 1:numel(poincareSection)
    poincareSection{idx}.numPoints = 80;  %#ok<SAGROW>
    % set maximum integration length to twice the expected circumference
    rOrbit = norm(poincareSection{idx}.endPosition(2,:)-poincareSection{idx}.endPosition(1,:));
    poincareSection{idx}.integrationLength = [0,2*(2*pi*rOrbit)]; %#ok<SAGROW>
end

% Plot Poincare sections
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),1:numel(poincareSection));
set(hPoincareSection,'color','g')
set(hPoincareSection,'linestyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'markerFaceColor','g')
drawnow

[shearline.etaPos,shearline.etaNeg] = lagrangian_shear(cgEigenvector,cgEigenvalue);
showGraph = true;
odeSolverOptions = odeset('relTol',1e-4);
dThresh = 1e-2;
closedOrbits = poincare_closed_orbit_multi(bickleyJet.flow,shearline,poincareSection,odeSolverOptions,dThresh,showGraph);
% Combine positive and negative eta closed orbits, and remove NaNs
closedOrbitsPos = cellfun(@(input)input{1},closedOrbits,'UniformOutput',false);
nanIdx = cellfun(@(input)any(isnan(input(:))),closedOrbitsPos);
closedOrbitsPos = closedOrbitsPos(~nanIdx);
closedOrbitsNeg = cellfun(@(input)input{2},closedOrbits,'UniformOutput',false);
nanIdx = cellfun(@(input)any(isnan(input(:))),closedOrbitsNeg);
closedOrbitsNeg = closedOrbitsNeg(~nanIdx);
closedOrbits = [closedOrbitsPos,closedOrbitsNeg];

hClosedOrbit = cellfun(@(input)plot(hAxes,input(:,1),input(:,2)),closedOrbits);
set(hClosedOrbit,'color','g')
drawnow

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,bickleyJet.strainline.maxLength,cgEigenvalue2,cgEigenvector1,bickleyJet.flow.domain,bickleyJet.flow.periodicBc);
for iClosedOrbit = 1:numel(closedOrbits)
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{iClosedOrbit});
end

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','r')
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor','r')

%% Attracting LCS analysis
bickleyJet.flow = set_flow_timespan([4*lengthX/u,0],bickleyJet.flow);

% Compute λ₂ and ξ₁
[cgEigenvalue,cgEigenvector] = eig_cgStrain(bickleyJet.flow,method,customEigMethod,coupledIntegration);
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(bickleyJet.flow.resolution));
cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(bickleyJet.flow.resolution),2]);

% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue2,diff(bickleyJet.flow.timespan));
hAxes = setup_figure(bickleyJet.flow.domain);
title(hAxes,'Attracting LCS')
plot_ftle(hAxes,bickleyJet.flow,ftle);
drawnow

% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection{1}.endPosition = [6.5,-1.4;4.5,-2]*1e6;
poincareSection{2}.endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection{3}.endPosition = [3.25,1.5;1.4,2.6]*1e6;
poincareSection{4}.endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection{5}.endPosition = [1.65e7,1.5e6;1.5e7,2.6e6];
for idx = 1:numel(poincareSection)
    poincareSection{idx}.numPoints = 80; %#ok<SAGROW>
    % set maximum integration length to twice the expected circumference
    rOrbit = norm(poincareSection{idx}.endPosition(2,:)-poincareSection{idx}.endPosition(1,:));
    poincareSection{idx}.integrationLength = [0,2*(2*pi*rOrbit)]; %#ok<SAGROW>
end

% Plot Poincare sections
hPoincareSection = arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1),poincareSection{idx}.endPosition(:,2)),1:numel(poincareSection));
set(hPoincareSection,'color','g')
set(hPoincareSection,'linestyle','--')
drawnow

[shearline.etaPos,shearline.etaNeg] = lagrangian_shear(cgEigenvector,cgEigenvalue);
showGraph = true;
odeSolverOptions = odeset('relTol',1e-4);
dThresh = 1e-2;
closedOrbits = poincare_closed_orbit_multi(bickleyJet.flow,shearline,poincareSection,odeSolverOptions,dThresh,showGraph);
% Combine positive and negative eta closed orbits, and remove NaNs
closedOrbitsPos = cellfun(@(input)input{1},closedOrbits,'UniformOutput',false);
nanIdx = cellfun(@(input)any(isnan(input(:))),closedOrbitsPos);
closedOrbitsPos = closedOrbitsPos(~nanIdx);
closedOrbitsNeg = cellfun(@(input)input{2},closedOrbits,'UniformOutput',false);
nanIdx = cellfun(@(input)any(isnan(input(:))),closedOrbitsNeg);
closedOrbitsNeg = closedOrbitsNeg(~nanIdx);
closedOrbits = [closedOrbitsPos,closedOrbitsNeg];

hClosedOrbit = cellfun(@(input)plot(hAxes,input(:,1),input(:,2)),closedOrbits);
set(hClosedOrbit,'color','g')
drawnow

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,bickleyJet.strainline.maxLength,cgEigenvalue2,cgEigenvector1,bickleyJet.flow.domain,bickleyJet.flow.periodicBc);
for iClosedOrbit = 1:numel(closedOrbits)
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{iClosedOrbit});
end

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color','r')
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor','r')
