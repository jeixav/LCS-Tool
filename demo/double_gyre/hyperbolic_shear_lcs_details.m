%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
flow.imposeIncompressibility = true;
flow.derivative = @(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega);
timespan = 20;

flow.domain = [-.1,2.1;-.05,1.05];
flow.timespan = [0,timespan];
flow.resolution = [551,276];

strainlineMaxLength = 20;

gridSpace = diff(flow.domain(1,:))/(double(flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;

forwardLcsColor = 'r';
backwardLcsColor = 'b';
shearLcsColor = [0,.6,0];

%% Forward-time LCS analysis
% Compute Cauchy-Green strain eigenvalues and eigenvectors
[flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow);
% FIXME Should use m-by-n or (m*n)-by-2 array forms throughout LCS Tool
cgEigenvalue = reshape(flow.cgEigenvalue,[fliplr(flow.resolution),2]);
cgEigenvector = reshape(flow.cgEigenvector,[fliplr(flow.resolution),4]);

% Plot finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(flow.timespan));
hAxes = setup_figure(flow.domain);
title(hAxes,'Forward-time LCS')
plot_ftle(hAxes,flow,ftle);
colormap(hAxes,flipud(gray))
drawnow

% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [.5,.6;.35,.5];
poincareSection(2).endPosition = [1.5,.4;1.7,.5];

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',shearLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',shearLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')
drawnow

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
closedOrbits = poincare_closed_orbit_multi(flow,shearline,poincareSection,'showGraph',true);

% Plot elliptic LCS
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbit
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)

% Plot all closed orbits of Poincare sections
for j = 1:nPoincareSection
    % η₊
    hOrbitsPos = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{1});
    set(hOrbitsPos,'color',shearLcsColor)
    % η₋
    hOrbitsNeg = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{2});
    set(hOrbitsNeg,'color',shearLcsColor)
end
drawnow

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,strainlineMaxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),flow.domain);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
end

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',forwardLcsColor)
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'MarkerSize',2)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor',forwardLcsColor)
drawnow

%% Backward-time LCS analysis
% Compute Cauchy-Green strain eigenvalues and eigenvectors
flow.timespan = [timespan,0];
[flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow);
cgEigenvalue = reshape(flow.cgEigenvalue,[fliplr(flow.resolution),2]);
cgEigenvector = reshape(flow.cgEigenvector,[fliplr(flow.resolution),4]);

% Plot backward time finite-time Lyapunov exponent
ftle = compute_ftle(cgEigenvalue(:,:,2),diff(flow.timespan));
hAxes = setup_figure(flow.domain);
title(hAxes,'Backward-time LCS')
plot_ftle(hAxes,flow,ftle);
colormap(hAxes,gray)
drawnow

% Define Poincare sections
poincareSection(1).endPosition = [.5,.6;.35,.5];
poincareSection(2).endPosition = [1.5,.4;1.7,.5];

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',shearLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',shearLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')
drawnow

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
closedOrbits = poincare_closed_orbit_multi(flow,shearline,poincareSection,'showGraph',true);

% Plot elliptic LCS
% η₊ outermost closed orbit
hClosedOrbitsEtaPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaPos,'color',shearLcsColor)
set(hClosedOrbitsEtaPos,'linewidth',2)
% η₋ outermost closed orbit
hClosedOrbitsEtaNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
set(hClosedOrbitsEtaNeg,'linewidth',2)

% Plot all closed orbits of Poincare sections
for j = 1:nPoincareSection
    % η₊
    hOrbitsPos = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{1});
    set(hOrbitsPos,'color',shearLcsColor)
    % η₋
    hOrbitsNeg = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{2});
    set(hOrbitsNeg,'color',shearLcsColor)
end
drawnow

% Compute strainlines
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,strainlineMaxLength,cgEigenvalue(:,:,2),cgEigenvector(:,:,1:2),flow.domain);

for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{1}{end});
    strainlinePosition = remove_strain_in_shear(strainlinePosition,closedOrbits{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
    idx = inpolygon(strainlineInitialPosition(1,:),strainlineInitialPosition(2,:),closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2));
    strainlineInitialPosition = strainlineInitialPosition(:,~idx);
end

% Plot strainlines
hStrainline = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainline,'color',backwardLcsColor)
hStrainlineInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainlineInitialPosition,'MarkerSize',2)
set(hStrainlineInitialPosition,'marker','o')
set(hStrainlineInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineInitialPosition,'MarkerFaceColor',backwardLcsColor)
