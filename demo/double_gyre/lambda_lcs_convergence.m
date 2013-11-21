function lambda_lcs_convergence

%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [-.1,2.1;-.05,1.05];
resolutionX = 650:50:900;
timespan = [0,20];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = true;

%% LCS parameters
% Lambda-lines
lambda = 1;
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [.5,.6;.35,.5];
poincareSection(2).endPosition = [1.5,.4;1.8,.5];

% Number of orbit seed points along each Poincare section
[poincareSection.numPoints] = deal(80);

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
    resolutionY = round(diff(domain(2,:))/gridSpace);
    resolution = [lResolutionX,resolutionY];
    disp(['Resolution: ',num2str(resolution)])
    
    hAxes = setup_figure(domain);
    title(hAxes,['Resolution: ',num2str(resolution(1)),'\times',num2str(resolution(2))])
    
    %% Cauchy-Green strain eigenvalues and eigenvectors
    [cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible);
    
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
    closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection);
    
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
