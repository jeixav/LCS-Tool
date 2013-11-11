%% Load data set
load('ocean_geostrophic_velocity.mat')

% Set velocity to zero at boundaries
vlon(:,[1,end],:) = 0;
vlon(:,:,[1,end]) = 0;
vlat(:,[1,end],:) = 0;
vlat(:,:,[1,end]) = 0;

strainlineLcsColor = 'r';
stretchlineLcsColor = 'b';
lambdaLineLcsColor = [0,.6,0];

%% Set parameters
% Define right hand side of ODE, ocean.flow.derivative
interpMethod = 'spline';
vlon_interpolant = griddedInterpolant({time,lat,lon},vlon,interpMethod);
vlat_interpolant = griddedInterpolant({time,lat,lon},vlat,interpMethod);
ocean.flow.derivative = @(t,y,~)flowdata_derivative(t,y,vlon_interpolant,vlat_interpolant);

% Set domain of initial conditions
% Center of domain [lon,lat]
center = [3,-31];
halfwidth = 3;
subdomain = [center(1)-halfwidth,center(1)+halfwidth;center(2)-halfwidth,center(2)+halfwidth];
ocean.flow = set_flow_domain(subdomain,ocean.flow);

% Set computation method for Cauchy-Green (CG) tensor
ocean.flow.cgStrainMethod.name = 'finiteDifference';
% Set if CG eigenvalues computed from main grid ('true' yields smoother eigenvalue fields)
ocean.flow.cgStrainMethod.eigenvalueFromMainGrid = false;
% Set auxiliary grid distance (relative value, i.e. 0.1 means 10% of maingrid size)
ocean.flow.cgStrainMethod.auxiliaryGridRelativeDelta = 0.1;
% Set computation method for eigenvectors
% false: use 'eig' function of MATLAB
% true: xi2 explicitly from auxiliary grid CG, xi1 as rotated xi2
ocean.flow.customEigMethod = false;
% Set if incompressibility of the flow is enforced,
% i.e., lambda1 = 1/lamda2
ocean.flow.imposeIncompressibility = true;
% Set resolution of subdomain
nxy = 400;
subdomainResolution = [nxy,nxy];
ocean.flow = set_flow_resolution(subdomainResolution,ocean.flow);
ocean.flow.timespan = [98,128];

lambda = 1;
shearlineOdeSolverOptions = odeset('relTol',1e-6);

strainlineOdeSolverOptions = odeset('relTol',1e-4);

gridSpace = diff(ocean.flow.domain(1,:))/(double(ocean.flow.resolution(1))-1);
localMaxDistance = 2*gridSpace;
hyperbolicLcsMaxLength = 20;

%% Cauchy-Green strain eigenvalues and eigenvectors
[ocean.flow.cgEigenvalue,ocean.flow.cgEigenvector] = eig_cgStrain(ocean.flow,ocean.flow.cgStrainMethod,ocean.flow.customEigMethod);

% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(ocean.flow.cgEigenvalue(:,2),fliplr(ocean.flow.resolution));
ftle_ = ftle(cgEigenvalue2,diff(ocean.flow.timespan));
hAxes = setup_figure(ocean.flow.domain);
title(hAxes,'Stretchline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')
plot_ftle(hAxes,ocean.flow,ftle_);
colormap(hAxes,flipud(gray))
drawnow

%% Shear LCSs
[ocean.shearline.etaPos,ocean.shearline.etaNeg] = lambda_line(ocean.flow.cgEigenvector,ocean.flow.cgEigenvalue,lambda);

% Define Poincare sections for closed orbit detection
% Poincare section should be placed with 1st point in center of elliptic region and
% with second point outside the elliptic region

poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});

% poincareSection(i).endPosition = [longitude1,latitude1;longitude2,latitude2]
poincareSection(1).endPosition = [3.15,-32.2;3.7,-31.6];
poincareSection(2).endPosition = [5,-31.6;5.3,-31.6];
poincareSection(3).endPosition = [4.8,-29.5;4.4,-29.5];
poincareSection(4).endPosition = [1.5,-30.9;1.9,-31.1];
poincareSection(5).endPosition = [2.9,-29.2;3.2,-29];
    
% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',lambdaLineLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')
drawnow

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length conservatively to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

% Closed orbit detection
disp('Detect elliptic LCS ...')
closedLambdaLine = poincare_closed_orbit_multi(ocean.flow,ocean.shearline,poincareSection,'odeSolverOptions',shearlineOdeSolverOptions,'showGraph',true);

% Plot lambda line LCSs
hOuterLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hOuterLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hOuterLambdaLineLcs = [hOuterLambdaLineLcsPos,hOuterLambdaLineLcsNeg];
set(hOuterLambdaLineLcs,'color',lambdaLineLcsColor)
set(hOuterLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaLinePos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{1});
    hClosedLambdaLineNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{2});
end
hClosedLambdaLine = horzcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
set(hClosedLambdaLine,'color',lambdaLineLcsColor)
drawnow

% Compute strainlines
disp('Detect hyperbolic LCS ...')
disp('Compute strainlines ...')
[strainlineLcs,strainlineLcsInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,hyperbolicLcsMaxLength,ocean.flow.cgEigenvalue(:,2),ocean.flow.cgEigenvector(:,1:2),ocean.flow.domain,ocean.flow.resolution,'odeSolverOptions',strainlineOdeSolverOptions);

% Plot hyperbolic strainline LCS
hStrainlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),strainlineLcs);
set(hStrainlineLcs,'color',strainlineLcsColor)
hStrainlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,strainlineLcsInitialPosition(1,idx),strainlineLcsInitialPosition(2,idx)),1:size(strainlineLcsInitialPosition,2));
set(hStrainlineLcsInitialPosition,'MarkerSize',5)
set(hStrainlineLcsInitialPosition,'marker','o')
set(hStrainlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStrainlineLcsInitialPosition,'MarkerFaceColor',strainlineLcsColor)

uistack(hPoincareSection,'top')
uistack(hLambdaLineLcs,'top')
uistack(hClosedLambdaLine,'top')
drawnow

%% Hyperbolic stretchline LCSs
% Plot finite-time Lyapunov exponent
hAxes = setup_figure(ocean.flow.domain);
title(hAxes,'Stretchline and \lambda-line LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')
plot_ftle(hAxes,ocean.flow,ftle_);
colormap(hAxes,flipud(gray))

% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
set(hPoincareSection,'color',lambdaLineLcsColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',lambdaLineLcsColor)
set(hPoincareSection,'MarkerEdgeColor','w')

% Plot lambda line LCSs
hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
set(hLambdaLineLcs,'color',lambdaLineLcsColor)
set(hLambdaLineLcs,'linewidth',2)

% Plot all closed lambda lines
hClosedLambdaLinePos = cell(nPoincareSection,1);
hClosedLambdaLineNeg = cell(nPoincareSection,1);
for j = 1:nPoincareSection
    hClosedLambdaLinePos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{1});
    hClosedLambdaLineNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{2});
end
hClosedLambdaLine = horzcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
set(hClosedLambdaLine,'color',lambdaLineLcsColor)
drawnow

stretchlineLcsOdeSolverOptions = odeset('relTol',1e-4);
[stretchlineLcs,stretchlineLcsInitialPosition] = seed_curves_from_lambda_max(localMaxDistance,hyperbolicLcsMaxLength,-ocean.flow.cgEigenvalue(:,1),ocean.flow.cgEigenvector(:,3:4),ocean.flow.domain,ocean.flow.resolution,'odeSolverOptions',stretchlineLcsOdeSolverOptions);

% Plot hyperbolic stretchline LCSs
hStretchlineLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchlineLcs);
set(hStretchlineLcs,'color',stretchlineLcsColor)
hStretchlineLcsInitialPosition = arrayfun(@(idx)plot(hAxes,stretchlineLcsInitialPosition(1,idx),stretchlineLcsInitialPosition(2,idx)),1:size(stretchlineLcsInitialPosition,2));
set(hStretchlineLcsInitialPosition,'MarkerSize',5)
set(hStretchlineLcsInitialPosition,'marker','o')
set(hStretchlineLcsInitialPosition,'MarkerEdgeColor','w')
set(hStretchlineLcsInitialPosition,'MarkerFaceColor',stretchlineLcsColor)

uistack(hPoincareSection,'top')
uistack(hLambdaLineLcs,'top')
uistack(hClosedLambdaLine,'top')
