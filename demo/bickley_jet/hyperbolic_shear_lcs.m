%% Input parameters
u = 62.66;
lengthX = pi*earthRadius;
lengthY = 1.77e6;
epsilon = [.075,.4,.3];

bickleyJet.flow.imposeIncompressibility = true;
bickleyJet.flow.periodicBc = [true,false];
perturbationCase = 3;
bickleyJet.flow = set_flow_derivative(@(t,x,useEoV)derivative(t,x,useEoV,u,lengthX,lengthY,epsilon,perturbationCase),bickleyJet.flow);
bickleyJet.flow = set_flow_resolution([500,200],bickleyJet.flow);
% magicNumber gives the domain an aspect ratio similar to that used in
% doi:10.1016/j.physd.2012.06.012 and ensures grid spacing is equal in the
% x and y directions. In doi:10.1016/j.physd.2012.06.012, 
% magicNumber = 2.2599.
magicNumber = .5*pi*earthRadius/lengthY*double(bickleyJet.flow.resolution(2)-1)/double(bickleyJet.flow.resolution(1)-1);
bickleyJet.flow = set_flow_domain([0,lengthX;[-1,1]*magicNumber*lengthY],bickleyJet.flow);
bickleyJet.flow = set_flow_timespan([0,4*lengthX/u],bickleyJet.flow);

shearlineOdeSolverOptions = odeset('relTol',1e-4);

bickleyJet.strainline = set_strainline_max_length(1e8);
gridSpace = diff(bickleyJet.flow.domain(1,:))/(double(bickleyJet.flow.resolution(1))-1);
localMaxDistance = 8*gridSpace;

forwardLcsColor = 'r';
backwardLcsColor = 'b';
shearLcsColor = [0,.6,0];

%% Forward-time LCS analysis
% Compute λ₂ and ξ₁
[cgEigenvalue,cgEigenvector] = eig_cgStrain(bickleyJet.flow);
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(bickleyJet.flow.resolution));
cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(bickleyJet.flow.resolution),2]);

% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [6.5,-1.4;5.2,-1.4]*1e6;
poincareSection(2).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection(3).endPosition = [3.25,1.5;1.4,2.6]*1e6;
poincareSection(4).endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection(5).endPosition = [1.65e7,1.5e6;1.5e7,2.6e6];

[poincareSection.numPoints] = deal(80);

nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

lambda = 1;
[shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedOrbits = poincare_closed_orbit_multi(bickleyJet.flow,shearline,poincareSection,'odeSolverOptions',shearlineOdeSolverOptions);

% Plot closed orbits
hAxes = setup_figure(bickleyJet.flow.domain);
title(hAxes,'Forward-time LCS')
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbits
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)
drawnow

% Compute strainlines
strainlinePosition = seed_curves_from_lambda_max(localMaxDistance,bickleyJet.strainline.maxLength,cgEigenvalue2,cgEigenvector1,bickleyJet.flow.domain,'periodicBc',bickleyJet.flow.periodicBc);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
end

% Plot hyperbolic LCS
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',forwardLcsColor)
drawnow

%% Backward-time LCS analysis
bickleyJet.flow = set_flow_timespan([4*lengthX/u,0],bickleyJet.flow);

% Compute λ₂ and ξ₁
[cgEigenvalue,cgEigenvector] = eig_cgStrain(bickleyJet.flow);
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(bickleyJet.flow.resolution));
cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(bickleyJet.flow.resolution),2]);

% Define Poincare sections
% Place first point in center of elliptic region and second point outside
% elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [6.5,-1.4;5.3,-1.4]*1e6;
poincareSection(2).endPosition = [1.35e7,-1.4e6;1.5e7,-.5e6];
poincareSection(3).endPosition = [3.25,1.5;1.4,2.6]*1e6;
poincareSection(4).endPosition = [1e7,1.5e6;8e6,2.6e6];
poincareSection(5).endPosition = [1.65e7,1.5e6;1.5e7,2.6e6];

[poincareSection.numPoints] = deal(80);

nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

[shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedOrbits = poincare_closed_orbit_multi(bickleyJet.flow,shearline,poincareSection,'odeSolverOptions',shearlineOdeSolverOptions);

% Plot closed orbits
hAxes = setup_figure(bickleyJet.flow.domain);
title(hAxes,'Backward-time LCS')
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbits
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)
drawnow

% Compute strainlines
strainlinePosition = seed_curves_from_lambda_max(localMaxDistance,bickleyJet.strainline.maxLength,cgEigenvalue2,cgEigenvector1,bickleyJet.flow.domain,'periodicBc',bickleyJet.flow.periodicBc);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
end

% Plot hyperbolic LCS
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',backwardLcsColor)
