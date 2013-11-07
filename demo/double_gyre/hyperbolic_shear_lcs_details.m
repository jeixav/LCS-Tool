%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;

flow.imposeIncompressibility = true;
flow.derivative = @(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega);
flow.domain = [-.1,2.1;-.05,1.05];
flow.timespan = [0,20];
flow.resolution = [551,276];

lambda = 1;

hyperbolicLcsMaxLength = 20;
gridSpace = diff(flow.domain(1,:))/(double(flow.resolution(1))-1);
strainlineLcsLocalMaxDistance = 2*gridSpace;

stretchlineLcsLocalMaxDistance = 4*gridSpace;

strainlineLcsColor = 'r';
stretchlineLcsColor = 'b';
lambdaLineLcsColor = [0,.6,0];

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvalue,cgEigenvector] = eig_cgStrain(flow);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(flow.resolution));
ftle_ = ftle(cgEigenvalue2,diff(flow.timespan));
hAxes = setup_figure(flow.domain);
title(hAxes,'Strainline and \lambda-line LCSs')
plot_ftle(hAxes,flow,ftle_);
colormap(hAxes,flipud(gray))
drawnow

%% Lambda-line LCSs
% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [.5,.6;.35,.5];
poincareSection(2).endPosition = [1.5,.4;1.7,.5];

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',neutralLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',neutralLcsColor)
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

[shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
closedOrbits = poincare_closed_orbit_multi(flow,shearline,poincareSection,'showGraph',true);

% Plot shear LCSs
% η₊ outermost closed orbit
hShearLcsPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hShearLcsPos,'color',lambdaLineLcsColor)
set(hShearLcsPos,'linewidth',2)
% η₋ outermost closed orbit
hShearLcsNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hShearLcsNeg,'color',lambdaLineLcsColor)
set(hShearLcsNeg,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaPos = cell(nPoincareSection,1);
hClosedLambdaNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaPos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{1});
    hClosedLambdaNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{2});
end
hClosedLambdaPos = horzcat(hClosedLambdaPos{:});
hClosedLambdaNeg = horzcat(hClosedLambdaNeg{:});
set(hClosedLambdaPos,'color',neutralLcsColor)
set(hClosedLambdaNeg,'color',neutralLcsColor)
drawnow

%% Hyperbolic strainline LCSs
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,hyperbolicLcsMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),flow.domain,flow.resolution);

% Plot hyperbolic strainline LCSs
hStrainLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainLcs,'color',unstableLcsColor)
hStrainLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainLcsInitialPosition,'MarkerSize',2)
set(hStrainLcsInitialPosition,'marker','o')
set(hStrainLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainLcsInitialPosition,'MarkerFaceColor',unstableLcsColor)

uistack(hShearLcsPos,'top')
uistack(hShearLcsNeg,'top')
uistack(hClosedLambdaPos,'top')
uistack(hClosedLambdaNeg,'top')
uistack(hPoincareSection,'top')
drawnow

%% Hyperbolic stretchline LCSs
% Plot finite-time Lyapunov exponent
hAxes = setup_figure(flow.domain);
title(hAxes,'Stretchline and \lambda-line LCSs')
plot_ftle(hAxes,flow,ftle_);
colormap(hAxes,flipud(gray))

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',neutralLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',neutralLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')

% Plot shear LCSs
% η₊ outermost closed orbit
hShearLcsPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
set(hShearLcsPos,'color',neutralLcsColor)
set(hShearLcsPos,'linewidth',2)
% η₋ outermost closed orbit
hShearLcsNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
set(hShearLcsNeg,'color',neutralLcsColor)
set(hShearLcsNeg,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaPos = cell(nPoincareSection,1);
hClosedLambdaNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaPos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{1});
    hClosedLambdaNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{2});
end
hClosedLambdaPos = horzcat(hClosedLambdaPos{:});
hClosedLambdaNeg = horzcat(hClosedLambdaNeg{:});
set(hClosedLambdaPos,'color',neutralLcsColor)
set(hClosedLambdaNeg,'color',neutralLcsColor)
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
[stretchlinePosition,stretchlineInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,hyperbolicLcsMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),flow.domain,flow.resolution);

% Plot hyperbolic stretchline LCSs
hStretchLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlinePosition);
set(hStretchLcs,'color',stableLcsColor)
hStretchLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineInitialPosition(1,idx),stretchlineInitialPosition(2,idx)),1:size(stretchlineInitialPosition,2));
set(hStretchLcsInitialPosition,'MarkerSize',2)
set(hStretchLcsInitialPosition,'marker','o')
set(hStretchLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchLcsInitialPosition,'MarkerFaceColor',stableLcsColor)

uistack(hShearLcsPos,'top')
uistack(hShearLcsNeg,'top')
uistack(hClosedLambdaPos,'top')
uistack(hClosedLambdaNeg,'top')
uistack(hPoincareSection,'top')
