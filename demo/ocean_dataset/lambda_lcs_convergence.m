function lambda_lcs_convergence

%% Input parameters
domain = [0,6;-34,-28];
resolutionX = 300:50:600;
timespan = [98,128];

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
cgAuxGridRelDelta = 0.1;

% Lambda-lines
lambda = 1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.15,-32.2;3.7,-31.6];
poincareSection(2).endPosition = [5,-31.6;5.3,-31.6];
poincareSection(3).endPosition = [4.8,-29.5;4.4,-29.5];
poincareSection(4).endPosition = [1.5,-30.9;1.9,-31.1];
poincareSection(5).endPosition = [2.9,-29.2;3.2,-29];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(100);

% Set maximum orbit length to twice the expected circumference
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
end

% Graphics properties
lambdaLineColor = [0,.6,0];

for m = 1:numel(resolutionX)
    % Make x and y grid spacing as equal as possible
    lResolutionX = resolutionX(m);
    gridSpace = diff(domain(1,:))/(double(lResolutionX)-1);
    resolutionY = round(diff(domain(2,:))/gridSpace) + 1;
    resolution = [lResolutionX,resolutionY];
    disp(['Resolution: ',num2str(resolution)])
    
    hAxes = setup_figure(domain);
    title(hAxes,['Resolution: ',num2str(resolution(1)),'\times',num2str(resolution(2))])
    xlabel(hAxes,'Longitude (\circ)')
    ylabel(hAxes,'Latitude (\circ)')

    %% Cauchy-Green strain eigenvalues and eigenvectors
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);
    
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
    
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
    closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions);
    
    % Plot lambda-line LCSs
    hLambdaLineLcsPos = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{1}{end}(:,1),closedLambdaLine{i}{1}{end}(:,2)),1:size(closedLambdaLine,2));
    hLambdaLineLcsNeg = arrayfun(@(i)plot(hAxes,closedLambdaLine{i}{2}{end}(:,1),closedLambdaLine{i}{2}{end}(:,2)),1:size(closedLambdaLine,2));
    hLambdaLineLcs = [hLambdaLineLcsPos,hLambdaLineLcsNeg];
    set(hLambdaLineLcs,'color',lambdaLineColor)
    set(hLambdaLineLcs,'linewidth',2)
    
    % Plot all closed lambda lines
    hClosedLambdaLinePos = cell(nPoincareSection,1);
    hClosedLambdaLineNeg = cell(nPoincareSection,1);
    for j = 1:nPoincareSection
        hClosedLambdaLinePos{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{1});
        hClosedLambdaLineNeg{j} = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),closedLambdaLine{j}{2});
    end
    hClosedLambdaLine = horzcat(hClosedLambdaLinePos{:},hClosedLambdaLineNeg{:});
    set(hClosedLambdaLine,'color',lambdaLineColor)
    uistack(hPoincareSection,'top')
    drawnow
end
