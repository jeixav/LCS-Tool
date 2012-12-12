% poincare_closed_orbit Find closed orbit using Poincare section
%
% DESCRIPTION
% This function is under development.
%
% EXAMPLE
% load('datasets/bickley_jet/bickleyJet1.mat')
% bickleyJet = shear_lcs_script(bickleyJet);
% set(findobj(gca,'tag','shearlineNegFiltered'),'visible','off')
% set(gca,'xlim',[3.4 4.8])
% set(gca,'ylim',[-1.3 -0.2])
% poincareSection.endPosition = [4.15 -.6; 4.15 -.3];
% poincareSection.numPoints = 30;
% odeSolverOptions = odeset('relTol',1e-6);
% poincare_closed_orbit(bickleyJet.flow,bickleyJet.shearline.etaPos,...
%   poincareSection,odeSolverOptions);

function closedOrbitPosition = poincare_closed_orbit(flow,vectorField,...
    poincareSection,odeSolverOptions)

% Initial positions for Poincare orbits
orbitInitialPositionX = linspace(poincareSection.endPosition(1,1),...
    poincareSection.endPosition(2,1),poincareSection.numPoints);
orbitInitialPositionY = linspace(poincareSection.endPosition(1,2),...
    poincareSection.endPosition(2,2),poincareSection.numPoints);
orbitInitialPosition = transpose([orbitInitialPositionX; ...
    orbitInitialPositionY]);

plot(orbitInitialPosition(:,1),orbitInitialPosition(:,2),'-x')

flowDomain = flow.domain;
flowResolution = flow.resolution;

% odeSolverOptions = odeset('RelTol',1e-6);
timespan = [0 5];
orbitPosition = cell(poincareSection.numPoints,1);

parfor idx = 1:poincareSection.numPoints
    orbitPosition{idx} = integrate_line_closed(timespan,...
        orbitInitialPosition(idx,:),flowDomain,flowResolution,...
        vectorField,odeSolverOptions);
end

arrayfun(@(idx)plot(orbitPosition{idx}(:,1),orbitPosition{idx}(:,2)),...
    1:poincareSection.numPoints)

orbitFinalPosition = cellfun(@(position)position(end,:),orbitPosition,...
    'UniformOutput',false);
orbitFinalPosition = cell2mat(orbitFinalPosition);

f1 = figure;
a1 = axes;
set(a1,'parent',f1);
set(a1,'nextplot','add')
set(a1,'box','on');
set(a1,'xgrid','on');
set(a1,'ygrid','on');
set(a1,'xlim',orbitInitialPosition([1,end],2))
plot(orbitInitialPosition(:,2),orbitFinalPosition(:,2) ...
    - orbitInitialPosition(:,2),'-x')
xlabel('s')
ylabel('p(s) - s')

closedOrbitPosition = [];
