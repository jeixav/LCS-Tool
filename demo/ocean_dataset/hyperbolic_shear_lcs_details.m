%% Input parameters
domain = [0,6;-34,-28];
resolution = [400,400];
timespan = [100,130];

%% Velocity definition
load('ocean_geostrophic_velocity.mat')
% Set velocity to zero at boundaries
vlon(:,[1,end],:) = 0;
vlon(:,:,[1,end]) = 0;
vlat(:,[1,end],:) = 0;
vlat(:,:,[1,end]) = 0;
interpMethod = 'spline';
vlon_interpolant = griddedInterpolant({time,lat,lon},vlon,interpMethod);
vlat_interpolant = griddedInterpolant({time,lat,lon},vlat,interpMethod);
lDerivative = @(t,x,~)flowdata_derivative(t,x,vlon_interpolant,vlat_interpolant);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgEigenvalueFromMainGrid = false;
cgAuxGridRelDelta = 0.01;

% Lambda-lines
lambdaLineOdeSolverOptions = odeset('relTol',1e-6,'initialStep',1e-2);
lambdaStep = 0.02;
lambdaRange = 0.90:lambdaStep:1.10;
% set flag to true to show Poincare return maps in closed orbit detection
showGraph = false;

% Strainlines
strainlineMaxLength = 20;
gridSpace = diff(domain(1,:))/(double(resolution(1))-1);
strainlineLocalMaxDistance = 2*gridSpace;
strainlineOdeSolverOptions = odeset('relTol',1e-6);

% Stretchlines
stretchlineMaxLength = 20;
stretchlineLocalMaxDistance = 4*gridSpace;
stretchlineOdeSolverOptions = odeset('relTol',1e-6);

% Graphics properties
strainlineColor = 'r';
stretchlineColor = 'b';
lambdaLineColor = [0,.6,0];
lcsInitialPositionMarkerSize = 2;

hAxes = setup_figure(domain);
title(hAxes,'Strainline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
ftle_(isnan(ftle_)) = max(ftle_(:));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
drawnow

save data.mat
% load data.mat

%% Lambda-line LCSs
% Define Poincare sections; first point in center of elliptic region and
% second point outside elliptic region
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

poincareSection(1).endPosition = [3.3,-32.1; 3.7,-31.6];
poincareSection(2).endPosition = [1.3,-30.9; 2.0,-31.2];
poincareSection(3).endPosition = [4.9,-29.6; 5.7,-29.6];
poincareSection(4).endPosition = [4.9,-31.4; 5.3,-31.4];
poincareSection(5).endPosition = [3.0,-29.3; 3.5,-29.3];

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',lambdaLineColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineColor)
set(hPoincareSection,'MarkerEdgeColor','w')
hNo = arrayfun(@(i)text(poincareSection(i).endPosition(1,1)-0.3,poincareSection(i).endPosition(1,2),num2str(i),'color','k','FontSize',16, 'fontname','lucida'),1:5);
drawnow

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    % FIXME Poincare section end point is not always close to vortex
    % centre, therefore rOrbit is not a good approximation of radius
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end

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
uistack(hNo,'top');
drawnow

%% Hyperbolic stretchline LCSs
hAxes = setup_figure(domain);
title(hAxes,'Stretchline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

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
hNo = arrayfun(@(i)text(poincareSection(i).endPosition(1,1)-0.3,poincareSection(i).endPosition(1,2),num2str(i),'color','k','FontSize',16, 'fontname','lucida'),1:5);

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
uistack(hNo,'top');

% print_eps(1,'LCS_strain');
% print_eps(2,'LCS_stretch');
