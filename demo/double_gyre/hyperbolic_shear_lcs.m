function hyperbolic_shear_lcs

%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;

flow.imposeIncompressibility = true;
flow.derivative = @(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega);
flow.domain = [-.1,2.1;-.05,1.05];
flow.timespan = [0,20];
flow.resolution = [551,276];

forwardLcsColor = 'r';
backwardLcsColor = 'b';
shearLcsColor = [0,.6,0];

%% Forward-time LCS analysis
% Compute Cauchy-Green strain eigenvalues and eigenvectors
[flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow);
% FIXME Should use m-by-n or (m*n)-by-2 array forms throughout LCS Tool
cgEigenvalue = reshape(flow.cgEigenvalue,[fliplr(flow.resolution),2]);
cgEigenvector = reshape(flow.cgEigenvector,[fliplr(flow.resolution),4]);

% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [.5,.6;.35,.5];
poincareSection(2).endPosition = [1.5,.4;1.7,.5];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(80);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

lambda = 1;
[shearline.etaPos,shearline.etaNeg] = lambda_line(flow.cgEigenvector,flow.cgEigenvalue,lambda);

% Compute closed shearlines
closedOrbits = poincare_closed_orbit_multi(flow,shearline,poincareSection);

% Plot elliptic LCS
hAxes = setup_figure(flow.domain);
title(hAxes,'Forward-time LCS')
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbit
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)
drawnow

% Compute strainlines
strainlineMaxLength = 20;
gridSpace = diff(flow.domain(1,:))/(double(flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;

strainlinePosition = seed_curves_from_lambda_max(localMaxDistance,strainlineMaxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),flow.domain);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
end

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',forwardLcsColor)
drawnow

%% Backward-time LCS analysis
% Compute Cauchy-Green strain eigenvalues and eigenvectors
flow.timespan = fliplr(flow.timespan);
[flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow);
cgEigenvalue = reshape(flow.cgEigenvalue,[fliplr(flow.resolution),2]);
cgEigenvector = reshape(flow.cgEigenvector,[fliplr(flow.resolution),4]);

% Define Poincare sections
poincareSection(1).endPosition = [.5,.6;.35,.5];
poincareSection(2).endPosition = [1.5,.4;1.7,.5];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(80);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

[shearline.etaPos,shearline.etaNeg] = lambda_line(flow.cgEigenvector,flow.cgEigenvalue,lambda);

% Compute closed shearlines
closedOrbits = poincare_closed_orbit_multi(flow,shearline,poincareSection);

% Plot elliptic LCS
hAxes = setup_figure(flow.domain);
title(hAxes,'Backward-time LCS')
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbit
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)
drawnow

% Compute strainlines
strainlinePosition = seed_curves_from_lambda_max(localMaxDistance,strainlineMaxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),flow.domain);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
end

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',backwardLcsColor)
