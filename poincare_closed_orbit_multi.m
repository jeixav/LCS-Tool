% poincare_closed_orbit_multi Find closed orbits of multiple Poincare
% sections
%
% SYNTAX
% [closedOrbits,orbits] = poincare_closed_orbit_multi(flow,shearline,...
%     PSList)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,
%     'odeSolverOptions',odeSolverOptions)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,
%     'dThresh',dThresh)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,
%     'showGraph',showGraph)
%
% INPUT ARGUMENTS
% PSList: 1-by-n struct of Poincare sections
% Format of PSList
% PSList(i).endPosition: [endPosition1x,endPosition1y;endPosition2x,...
%     endPosition2y];
% PSList(i).numPoints: number of initial positions along Poincare section
% from which closed orbit candidates will be launched
% PSList(i).orbitMaxLength: maximum length allowed for closed orbits to
% limit integration time.
% dThresh: poincare_closed_orbit dThresh input; default is 1e-2.
% showGraph: logical value to control display of Poincare section return
% map plots; default is false.
%
% OUTPUT
% closedOrbits{}{}: Positions of closed orbits
% Format of closeOrbits
% closedOrbits{i}{1}{1}: innermost closed orbit around Poincare section i
% in etaPos field
% closedOrbits{i}{2}{end}: outermost closed orbit around Poincare section i
% in etaNeg field
% orbits{}{}{}: Positions of all orbits
% Format: orbits{1}{2}{3}: 3rd {3} orbit of 1st {1} Poincare section in
% etaNeg {2} field

function [closedOrbits,orbits] = poincare_closed_orbit_multi(flow,shearline,PSList,varargin)

p = inputParser;
addRequired(p,'flow',@isstruct)
addRequired(p,'shearline',@isstruct)
addRequired(p,'PSList',@isstruct)
addParamValue(p,'dThresh',1e-2,@(dThresh)validateattributes(dThresh,{'double'},{'scalar'}))
addParamValue(p,'odeSolverOptions',odeset)
addParamValue(p,'showGraph',false,@islogical)
parse(p,flow,shearline,PSList,varargin{:})
dThresh = p.Results.dThresh;
odeSolverOptions = p.Results.odeSolverOptions;
showGraph = p.Results.showGraph;

nPoincareSection = numel(PSList);
closedOrbits = cell(1,nPoincareSection);
orbits = cell(1,nPoincareSection);

nBisection = 5;

for i = 1:nPoincareSection
    % define current Poincare section
    poincareSection.endPosition = PSList(i).endPosition;
    poincareSection.numPoints = PSList(i).numPoints;
    poincareSection.integrationLength = [0,PSList(i).orbitMaxLength];
    
    [closedOrbitsPos,orbitsPos] = poincare_closed_orbit(flow,shearline.etaPos,poincareSection,odeSolverOptions,nBisection,dThresh,showGraph);
    closedOrbits{i}{1} = closedOrbitsPos;
    orbits{i}{1} = orbitsPos;
    
    [closedOrbitsNeg,orbitsNeg] = poincare_closed_orbit(flow,shearline.etaNeg,poincareSection,odeSolverOptions,nBisection,dThresh,showGraph);
    closedOrbits{i}{2} = closedOrbitsNeg;
    orbits{i}{2} = orbitsNeg;    
end
