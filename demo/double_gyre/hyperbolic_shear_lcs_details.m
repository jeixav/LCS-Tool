%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolution = [500,250];
timespan = [0,10];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% Lambda-lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [0.55,0.55;0.2,0.2];
poincareSection(2).endPosition = [1.53,.45;1.95,0.1];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
lambdaStep = 0.01;
lambdaRange = 0.95:lambdaStep:1.05;
% set flag to true to show Poincare return maps in closed orbit detection
showGraph = false;

% Strainlines
strainlineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLocalMaxDistance = 2*gridSpace;
strainlineOdeSolverOptions = odeset('relTol',1e-6);

% Stretchlines
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 10*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-6);

% Graphics properties
strainlineColor = 'r';
stretchlineColor = 'b';
lambdaLineColor = [0,.6,0];
lcsInitialPositionMarkerSize = 2;

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs');

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
drawnow

%% Lambda-line LCSs
% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',lambdaLineColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineColor)
set(hPoincareSection,'MarkerEdgeColor','w')
drawnow

% find closed orbits for range of lambda values
closedLambdaLineArea = zeros(1,nPoincareSection);
lambda0 = nan(1,nPoincareSection);
k=0;
for lambda = lambdaRange
    k=k+1;
    
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
    shearline.etaPos = real(shearline.etaPos);
    shearline.etaNeg = real(shearline.etaNeg);      
    
    closedLambdaLineCandidate = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',showGraph);
    
    % keep outermost closed orbit
    for i = 1:nPoincareSection
        for j=1:2 % etaPos,etaNeg
            orbitArea(j) = polyarea(closedLambdaLineCandidate{i}{j}{end}(:,1),closedLambdaLineCandidate{i}{j}{end}(:,2));
        end        
        if max(orbitArea) > closedLambdaLineArea(i)
            closedLambdaLineArea(i) = max(orbitArea);
            closedLambdaLine{i}{1}{1} = closedLambdaLineCandidate{i}{1}{1};
            closedLambdaLine{i}{2}{1} = closedLambdaLineCandidate{i}{2}{1};
            % keep lambda values associated to closed orbits            
            lambda0(i) = lambda;
        end        
    end    
end

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineColor)
set(hLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for i = 1:nPoincareSection
    hClosedLambdaLinePos{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{1});
    hClosedLambdaLineNeg{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{2});
end
hClosedLambdaLine = vertcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
set(hClosedLambdaLine,'color',lambdaLineColor)
drawnow

%% Hyperbolic strainline LCSs
[strainlineLcs,strainlineLcsInitialPosition] = seed_curves_from_lambda_max(strainlineLocalMaxDistance,strainlineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',strainlineOdeSolverOptions);

% Remove strainlines inside elliptic regions
for i = 1:nPoincareSection
    % Remove strainlines inside elliptic regions
    strainlineLcs = remove_strain_in_shear(strainlineLcs,closedLambdaLine{i}{1}{end});
    strainlineLcs = remove_strain_in_shear(strainlineLcs,closedLambdaLine{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(strainlineLcsInitialPosition(1,:),strainlineLcsInitialPosition(2,:),closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2));
    strainlineLcsInitialPosition = strainlineLcsInitialPosition(:,~idx);
    idx = inpolygon(strainlineLcsInitialPosition(1,:),strainlineLcsInitialPosition(2,:),closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2));
    strainlineLcsInitialPosition = strainlineLcsInitialPosition(:,~idx);
end

% Plot hyperbolic strainline LCSs
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs);
set(hStrainlineLcs,'color',strainlineColor)
hStrainlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineLcsInitialPosition(1,idx),strainlineLcsInitialPosition(2,idx)),1:size(strainlineLcsInitialPosition,2));
set(hStrainlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStrainlineLcsInitialPosition,'marker','o')
set(hStrainlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineLcsInitialPosition,'MarkerFaceColor',strainlineColor)
uistack(hLambdaLineLcs,'top')
uistack(hClosedLambdaLine,'top')
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
set(hPoincareSection,'color',lambdaLineColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineColor)
set(hPoincareSection,'MarkerEdgeColor','w')

% Plot lambda-line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineColor)
set(hLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for i = 1:nPoincareSection
    hClosedLambdaLinePos{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{1});
    hClosedLambdaLineNeg{i} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{i}{2});
end
hClosedLambdaLine = vertcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
set(hClosedLambdaLine,'color',lambdaLineColor)
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minimums
[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(stretchlineLocalMaxDistance,stretchlineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchlineOdeSolverOptions);

% Remove stretchlines inside elliptic regions
for i = 1:nPoincareSection
    % Remove stretchlines inside elliptic regions
    stretchlineLcs = remove_strain_in_shear(stretchlineLcs,closedLambdaLine{i}{1}{end});
    stretchlineLcs = remove_strain_in_shear(stretchlineLcs,closedLambdaLine{i}{2}{end});
    % Remove initial positions inside elliptic regions
    idx = inpolygon(stretchlineLcsInitialPosition(1,:),stretchlineLcsInitialPosition(2,:),closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2));
    stretchlineLcsInitialPosition = stretchlineLcsInitialPosition(:,~idx);
    idx = inpolygon(stretchlineLcsInitialPosition(1,:),stretchlineLcsInitialPosition(2,:),closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2));
    stretchlineLcsInitialPosition = stretchlineLcsInitialPosition(:,~idx);
end

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs);
set(hStretchlineLcs,'color',stretchlineColor)
hStretchlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineLcsInitialPosition(1,idx),stretchlineLcsInitialPosition(2,idx)),1:size(stretchlineLcsInitialPosition,2));
set(hStretchlineLcsInitialPosition,'MarkerSize',lcsInitialPositionMarkerSize)
set(hStretchlineLcsInitialPosition,'marker','o')
set(hStretchlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchlineLcsInitialPosition,'MarkerFaceColor',stretchlineColor)
uistack(hLambdaLineLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')

print_eps(1,'LCS_strain');
print_eps(2,'LCS_stretch');
