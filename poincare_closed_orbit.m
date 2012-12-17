% poincare_closed_orbit Find closed orbit using Poincare section
%
% DESCRIPTION
% This function is under development.
%
% EXAMPLE
% addpath('flow_templates')
% matlabpool('open')
% pctRunOnAll javaaddpath('ParforProgress2')
% 
% bickleyJet = bickley_jet(1);
% 
% bickleyJet.flow = set_flow_resolution(40,bickleyJet.flow);
% bickleyJet.flow = set_flow_ode_solver_options(odeset('relTol',1e-4,...
%     'absTol',1e-6),bickleyJet.flow);
% bickleyJet.flow = set_flow_domain([3.5 4.7; -1.3 -.2],bickleyJet.flow);
% bickleyJet.flow.imposeIncompressibility = false;
% 
% bickleyJet.shearline = rmfield(bickleyJet.shearline,'odeSolverOptions');
% bickleyJet.shearline = set_shearline_resolution([4 4],...
%     bickleyJet.shearline);
% 
% showPlot.shearlinePosFiltered = false;
% showPlot.shearlineNegFiltered = false;
% showPlot.etaPosQuiver = true;
% 
% bickleyJet = shear_lcs_script(bickleyJet,showPlot);
% 
% poincareSection.endPosition = [4.15 -.6; 4.15 -.3];
% poincareSection.numPoints = 50;
% odeSolverOptions = odeset('relTol',1e-6);
% 
% closedOrbitInitialPositionY = poincare_closed_orbit(bickleyJet.flow,...
%     bickleyJet.shearline.etaPos,poincareSection,odeSolverOptions);
% 
% disp('Closed orbit positions:')
% disp(num2str(transpose(closedOrbitInitialPositionY)))

function closedOrbitInitialPositionY = poincare_closed_orbit(flow,...
    vectorField,poincareSection,odeSolverOptions,showGraph)

narginchk(4,5)

if nargin == 4
    showGraph = true;
end

% Initial positions for Poincare orbits
orbitInitialPositionX = linspace(poincareSection.endPosition(1,1),...
    poincareSection.endPosition(2,1),poincareSection.numPoints);
orbitInitialPositionY = linspace(poincareSection.endPosition(1,2),...
    poincareSection.endPosition(2,2),poincareSection.numPoints);
orbitInitialPosition = transpose([orbitInitialPositionX; ...
    orbitInitialPositionY]);

if showGraph
    hPoincareSection = plot(orbitInitialPosition(:,1),...
        orbitInitialPosition(:,2),'-x');
    hParent = get(hPoincareSection,'parent');
end

flowDomain = flow.domain;
flowResolution = flow.resolution;

% FIXME The timespan should not be hard-coded
timespan = [0 5];
orbitPosition = cell(poincareSection.numPoints,1);

parfor idx = 1:poincareSection.numPoints
    orbitPosition{idx} = integrate_line_closed(timespan,...
        orbitInitialPosition(idx,:),flowDomain,flowResolution,...
        vectorField,poincareSection.endPosition,odeSolverOptions);
end

if showGraph
    hOrbit = arrayfun(@(idx)plot(orbitPosition{idx}(:,1),...
        orbitPosition{idx}(:,2)),1:poincareSection.numPoints);
end

orbitFinalPosition = cellfun(@(position)position(end,:),orbitPosition,...
    'UniformOutput',false);
orbitFinalPosition = cell2mat(orbitFinalPosition);

if showGraph
    hfigure = figure;
    hAxes = axes;
    set(hAxes,'parent',hfigure);
    set(hAxes,'nextplot','add')
    set(hAxes,'box','on');
    set(hAxes,'xgrid','on');
    set(hAxes,'ygrid','on');
    set(hAxes,'xlim',orbitInitialPosition([1,end],2))
    plot(hAxes,orbitInitialPosition(:,2),orbitFinalPosition(:,2) ...
        - orbitInitialPosition(:,2),'-x')
    xlabel(hAxes,'s')
    ylabel(hAxes,'p(s) - s')
end

[~,closedOrbitInitialPositionY] = crossing(orbitFinalPosition(:,2) ...
    - orbitInitialPosition(:,2),orbitInitialPosition(:,2));

if showGraph
    nClosedOrbit = numel(closedOrbitInitialPositionY);
    plot(hAxes,closedOrbitInitialPositionY,...
        zeros(1,nClosedOrbit),'ro')

    % FIXME Vertical poincare section assumed
    closedOrbitInitialPositionX = poincareSection.endPosition(1,1).*...
        ones(1,nClosedOrbit);
    orbitInitialPosition = transpose([closedOrbitInitialPositionX; ...
        closedOrbitInitialPositionY]);
    parfor idx = 1:nClosedOrbit
        closedOrbitPosition{idx} = integrate_line_closed(timespan,...
            orbitInitialPosition(idx,:),flowDomain,flowResolution,...
            vectorField,poincareSection.endPosition,odeSolverOptions);
    end

    hClosedOrbit = arrayfun(@(idx)plot(hParent,...
        closedOrbitPosition{idx}(:,1),closedOrbitPosition{idx}(:,2)),...
        1:nClosedOrbit);
    set(hClosedOrbit,'color','black')
    set(hClosedOrbit,'linewidth',2)
    
end
