%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;

incompressible = true;
derivative = @(t,x,useEoV)derivative(t,x,useEoV,epsilon,amplitude,omega);
domain = [-.1,2.1;-.05,1.05];
timespan = [0,20];
resolution = [551,276];

lambda = 1;

hyperbolicLcsMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLcsLocalMaxDistance = 2*gridSpace;

stretchlineLcsLocalMaxDistance = 4*gridSpace;

strainlineLcsColor = 'r';
stretchlineLcsColor = 'b';
lambdaLineLcsColor = [0,.6,0];

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvalue,cgEigenvector] = eig_cgStrain(derivative,domain,timespan,resolution,'incompressible',incompressible);

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
closedOrbits = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'showGraph',true);

% Plot lambda-line LCSs
% η₊ outermost closed lambda-lines
hShearLcsPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
% η₋ outermost closed lambda-lines
hShearLcsNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
hShearLcs = [hShearLcsPos,hShearLcsNeg];
set(hShearLcs,'color',lambdaLineLcsColor)
set(hShearLcs,'linewidth',2)

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

%% Hyperbolic strainline LCSs
[strainlinePosition,strainlineInitialPosition] = seed_curves_from_lambda_max(strainlineLcsLocalMaxDistance,hyperbolicLcsMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution);

% Plot hyperbolic strainline LCSs
hStrainLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlinePosition);
set(hStrainLcs,'color',strainlineLcsColor)
hStrainLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineInitialPosition(1,idx),strainlineInitialPosition(2,idx)),1:size(strainlineInitialPosition,2));
set(hStrainLcsInitialPosition,'MarkerSize',2)
set(hStrainLcsInitialPosition,'marker','o')
set(hStrainLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainLcsInitialPosition,'MarkerFaceColor',strainlineLcsColor)

uistack(hShearLcs,'top')
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
% η₊ outermost closed lambda-lines
hShearLcsPos = arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
% η₋ outermost closed lambda-lines
hShearLcsNeg = arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
hShearLcs = [hShearLcsPos,hShearLcsNeg];
set(hShearLcs,'color',lambdaLineLcsColor)
set(hShearLcs,'linewidth',2)

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
[stretchlinePosition,stretchlineInitialPosition] = seed_curves_from_lambda_max(stretchlineLcsLocalMaxDistance,hyperbolicLcsMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution);

% Plot hyperbolic stretchline LCSs
hStretchLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlinePosition);
set(hStretchLcs,'color',stretchlineLcsColor)
hStretchLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineInitialPosition(1,idx),stretchlineInitialPosition(2,idx)),1:size(stretchlineInitialPosition,2));
set(hStretchLcsInitialPosition,'MarkerSize',2)
set(hStretchLcsInitialPosition,'marker','o')
set(hStretchLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchLcsInitialPosition,'MarkerFaceColor',stretchlineLcsColor)

uistack(hShearLcs,'top')
uistack(hClosedLambda,'top')
uistack(hPoincareSection,'top')
