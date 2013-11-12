%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [-.1,2.1;-.05,1.05];
resolution = [551,276];
timespan = [0,20];

%% Velocity definition
derivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
% Lambda-lines
lambda = 1;

% Strainlines
strainlineLcsMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLcsLocalMaxDistance = 2*gridSpace;

% Stretchlines
stretchlineLcsMaxLength = 20;
stretchlineLcsLocalMaxDistance = 4*gridSpace;

% Graphics properties
strainlineLcsColor = 'r';
stretchlineLcsColor = 'b';
lambdaLineLcsColor = [0,.6,0];

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvalue,cgEigenvector] = eig_cgStrain(derivative,domain,resolution,timespan,'incompressible',incompressible);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
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
set(hPoincareSection,'color',lambdaLineLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineLcsColor)
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
closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'showGraph',true);

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineLcsColor)
set(hLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaPos = cell(nPoincareSection,1);
hClosedLambdaNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaPos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{1});
    hClosedLambdaNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{2});
end
hClosedLambda = horzcat(hClosedLambdaPos{:},hClosedLambdaNeg{:});
set(hClosedLambda,'color',lambdaLineLcsColor)
drawnow

%% Hyperbolic strainline LCSs
[strainlineLcs,strainlineInitialPosition] = seed_curves_from_lambda_max(strainlineLcsLocalMaxDistance,strainlineLcsMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution);

% Plot hyperbolic strainline LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs);
set(hStrainlineLcs,'color',strainlineLcsColor)
hStrainLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainLcsInitialPosition,'MarkerSize',2)
set(hStrainLcsInitialPosition,'marker','o')
set(hStrainLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainLcsInitialPosition,'MarkerFaceColor',strainlineLcsColor)

uistack(hLambdaLineLcs,'top')
uistack(hClosedLambda,'top')
uistack(hPoincareSection,'top')
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(domain);
title(hAxes,'Stretchline and \lambda-line LCSs')

% Plot finite-time Lyapunov exponent
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',lambdaLineLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineLcsColor)
set(hLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaPos = cell(nPoincareSection,1);
hClosedLambdaNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaPos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{1});
    hClosedLambdaNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedOrbits{j}{2});
end
hClosedLambda = horzcat(hClosedLambdaPos{:},hClosedLambdaNeg{:});
set(hClosedLambda,'color',lambdaLineLcsColor)
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
[stretchlineLcs,stretchlineInitialPosition] = seed_curves_from_lambda_max(stretchlineLcsLocalMaxDistance,stretchlineLcsMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution);

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs);
set(hStretchlineLcs,'color',stretchlineLcsColor)
hStretchLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineInitialPosition(1,idx),stretchlineInitialPosition(2,idx)),1:size(stretchlineInitialPosition,2));
set(hStretchLcsInitialPosition,'MarkerSize',2)
set(hStretchLcsInitialPosition,'marker','o')
set(hStretchLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchLcsInitialPosition,'MarkerFaceColor',stretchlineLcsColor)

uistack(hLambdaLineLcs,'top')
uistack(hClosedLambda,'top')
uistack(hPoincareSection,'top')
