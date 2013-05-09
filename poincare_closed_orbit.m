% poincare_closed_orbit Find closed orbit using Poincare section
%
% DESCRIPTION
% This function is under development.
%
% SYNTAX
% closedOrbitInitialPosition = poincare_closed_orbit(flow,vectorField,poincareSection,odeSolverOptions,showGraph)
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
% [bickleyJet,hAxes] = shear_lcs_script(bickleyJet,showPlot);
% 
% poincareSection.endPosition = [4.15,-.6;4.15,-.3];
% poincareSection.numPoints = 50;
% odeSolverOptions = odeset('relTol',1e-6);
% 
% closedOrbitInitialPosition = poincare_closed_orbit(bickleyJet.flow,...
%     bickleyJet.shearline.etaPos,poincareSection,odeSolverOptions,hAxes);
% 
% disp('Closed orbit positions:')
% disp(num2str(transpose(closedOrbitInitialPosition)))

function [closedOrbitInitialPosition,closedOrbitPosition,hAxes] = poincare_closed_orbit(flow,vectorField,poincareSection,odeSolverOptions,showGraph)

narginchk(4,5)

if nargin == 4
    showGraph = false;
else
    hAxes = showGraph;
end

% Initial positions for Poincare orbits
orbitInitialPositionX = linspace(poincareSection.endPosition(1,1),poincareSection.endPosition(2,1),poincareSection.numPoints);
orbitInitialPositionY = linspace(poincareSection.endPosition(1,2),poincareSection.endPosition(2,2),poincareSection.numPoints);
orbitInitialPosition = transpose([orbitInitialPositionX;orbitInitialPositionY]);

if showGraph
%     hPoincareSection = plot(orbitInitialPosition(:,1),...
%         orbitInitialPosition(:,2),'-x');
%     set(hPoincareSection,'parent',hAxes)
    hParent = hAxes;
end

flowDomain = flow.domain;
flowResolution = flow.resolution;

timespan = [0,poincareSection.tmax];
orbitPosition = cell(poincareSection.numPoints,1);

parfor idx = 1:poincareSection.numPoints
    orbitPosition{idx} = integrate_line_closed(timespan,orbitInitialPosition(idx,:),flowDomain,flowResolution,vectorField,poincareSection.endPosition,odeSolverOptions);
end

% if showGraph
%     arrayfun(@(idx)plot(hParent,orbitPosition{idx}(:,1),...
%         orbitPosition{idx}(:,2)),1:poincareSection.numPoints);
% end

orbitFinalPosition = cellfun(@(position)position(end,:),orbitPosition,'UniformOutput',false);
orbitFinalPosition = cell2mat(orbitFinalPosition);

if showGraph
    hfigure = figure;
    hAxes = axes('parent',hfigure);
    set(hAxes,'nextplot','add')
    set(hAxes,'box','on');
    set(hAxes,'xgrid','on');
    set(hAxes,'ygrid','on');
    
    xLength = sqrt(diff(poincareSection.endPosition(:,1))^2 + diff(poincareSection.endPosition(:,2))^2);
    theta = atan((poincareSection.endPosition(1,2) - poincareSection.endPosition(2,2))/(poincareSection.endPosition(1,1) - poincareSection.endPosition(2,1)));
    rotationMatrix = [cos(theta),-sin(theta);sin(theta),cos(theta)];

    % Translate to origin
    s(:,1) = orbitInitialPosition(:,1) - poincareSection.endPosition(1,1);
    s(:,2) = orbitInitialPosition(:,2) - poincareSection.endPosition(1,2);

    t(:,1) = orbitFinalPosition(:,1) - poincareSection.endPosition(1,1);
    t(:,2) = orbitFinalPosition(:,2) - poincareSection.endPosition(1,2);

    % Rotate to 0
    s = rotationMatrix*transpose(s);
    t = rotationMatrix*transpose(t);
    
    s = transpose(s);
    t = transpose(t);
    
    set(hAxes,'xlim',[0 xLength])
    plot(hAxes,s,t-s,'k')
    xlabel(hAxes,'s')
    ylabel(hAxes,'p(s) - s')
end

[~,closedOrbitInitialPosition] = crossing(t(:,1) - s(:,1),s(:,1));

if showGraph
    nClosedOrbit = numel(closedOrbitInitialPosition);
    hPlot = plot(hAxes,closedOrbitInitialPosition,zeros(1,nClosedOrbit));
    set(hPlot,'linestyle','none')
    set(hPlot,'marker','o')
    set(hPlot,'MarkerFaceColor','b')
    set(hPlot,'MarkerEdgeColor','b')
    set(hPlot,'tag','PoincareZero')
    
    % Rotate to theta
    xx = [transpose(closedOrbitInitialPosition),zeros(numel(closedOrbitInitialPosition),1)];
    xx = rotationMatrix\transpose(xx);
    xx = transpose(xx);
    
    % Translate from origin
    closedOrbitInitialPositionX = xx(:,1) + poincareSection.endPosition(1,1);
    closedOrbitInitialPositionY = xx(:,2) + poincareSection.endPosition(1,2);
    
    closedOrbitInitialPosition = [closedOrbitInitialPositionX,closedOrbitInitialPositionY];
    parfor idx = 1:nClosedOrbit
        closedOrbitPosition{idx} = integrate_line_closed(timespan,closedOrbitInitialPosition(idx,:),flowDomain,flowResolution,vectorField,poincareSection.endPosition,odeSolverOptions);
    end
    
    if ~isempty(closedOrbitInitialPositionY)
        hClosedOrbit = arrayfun(@(idx)plot(hParent,closedOrbitPosition{idx}(:,1),closedOrbitPosition{idx}(:,2)),1:nClosedOrbit);
        set(hClosedOrbit,'color','black')
        set(hClosedOrbit,'linewidth',2)
        set(hClosedOrbit,'tag','closedOrbit')
    end
    
end
