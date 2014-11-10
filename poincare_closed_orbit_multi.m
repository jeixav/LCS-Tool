% poincare_closed_orbit_multi Find closed orbits of multiple Poincare
% sections
%
% SYNTAX
% [closedOrbits,orbits] = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,PSList)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,'nBisection',nBisection)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,'dThresh',dThresh)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,'odeSolverOptions',odeSolverOptions)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,'periodicBc',periodicBc)
% [closedOrbits,orbits] = poincare_closed_orbit_multi(...,'showGraph',showGraph)
% [closedOrbits,orbits,hFigure] = poincare_closed_orbit_multi(...,'showGraph',showGraph)
%
% INPUT ARGUMENTS
% PSList: 1-by-n struct of Poincare sections
%
% PSList(i).endPosition: [endPosition1x,endPosition1y;endPosition2x,...
% endPosition2y];
%
% PSList(i).numPoints: number of initial positions along Poincare section
% from which closed orbit candidates will be launched
%
% PSList(i).orbitMaxLength: maximum length allowed for closed orbits to
% limit integration time.
%
% nBisection: poincare_closed_orbit nBisection input; default is 5.
%
% dThresh: poincare_closed_orbit dThresh input; default is 1e-2.
%
% showGraph: logical value to control display of Poincare section return
% map plots; default is false.
%
% OUTPUT ARGUMENTS
% closedOrbits: Closed orbit positions; cell array of size(PSList)
% closedOrbits{m}{1}{1}: innermost closed orbit around Poincare section m
% in etaPos field
%
% closedOrbits{m}{1}{end}: outermost closed orbit around Poincare section m
% in etaPos field
%
% closedOrbits{m}{2}{end}: outermost closed orbit around Poincare section m
% in etaNeg field
%
% orbits: Positions of all orbits
% orbits{1}{2}{3}: 3rd {3} orbit of 1st {1} Poincare section in etaNeg {2}
% field

function [closedOrbits,orbits,varargout] = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,PSList,varargin)

nargoutchk(1,3)

p = inputParser;

% FIXME Make validationFcn common with eig_cgStrain and poincare_closed_orbit_range
addRequired(p,'domain',@(input)validateattributes(input,{'double'},{'size',[2,2],'real','finite'}))
addRequired(p,'resolution',@(input)validateattributes(input,{'double'},{'size',[1,2],'real','finite'}))
addRequired(p,'etaPos',@(input)validateattributes(input,{'double'},{'size',[prod(resolution),2]}))
addRequired(p,'etaNeg',@(input)validateattributes(input,{'double'},{'size',[prod(resolution),2]}))
addRequired(p,'PSList',@isstruct)
addParameter(p,'nBisection',5,@(input)validateattributes(input,{'numeric'},{'scalar','integer','positive'}))
addParameter(p,'dThresh',1e-2,@(input)validateattributes(input,{'double'},{'scalar'}))
addParameter(p,'odeSolverOptions',odeset)
addParameter(p,'periodicBc',[false,false],@(input)validateattributes(input,{'logical'},{'size',[1,2]}));
addParameter(p,'showGraph',false,@islogical)

parse(p,domain,resolution,etaPos,etaNeg,PSList,varargin{:})

nBisection = p.Results.nBisection;
dThresh = p.Results.dThresh;
odeSolverOptions = p.Results.odeSolverOptions;
periodicBc = p.Results.periodicBc;
showGraph = p.Results.showGraph;

nPoincareSection = numel(PSList);
closedOrbits = cell(1,nPoincareSection);
orbits = cell(1,nPoincareSection);

if showGraph
    hFigure = nan(nPoincareSection,2);
end

for i = 1:nPoincareSection
    % define current Poincare section
    poincareSection.endPosition = PSList(i).endPosition;
    poincareSection.numPoints = PSList(i).numPoints;
    poincareSection.integrationLength = [0,PSList(i).orbitMaxLength];
    
    if showGraph
        [closedOrbitsPos,orbitsPos,hFigure(i,1)] = poincare_closed_orbit(domain,resolution,etaPos,poincareSection,'odeSolverOptions',odeSolverOptions,'dThresh',dThresh,'nBisection',nBisection,'periodicBc',periodicBc,'showGraph',showGraph);
        hTitle = get(get(hFigure(i,1),'CurrentAxes'),'Title');
        titleString = get(hTitle,'String');
        newTitleString = sprintf('%s #%u \\eta_+',titleString,i);
        set(hTitle,'String',newTitleString)
    else
        [closedOrbitsPos,orbitsPos] = poincare_closed_orbit(domain,resolution,etaPos,poincareSection,'odeSolverOptions',odeSolverOptions,'dThresh',dThresh,'nBisection',nBisection,'periodicBc',periodicBc,'showGraph',showGraph);
    end
    closedOrbits{i}{1} = closedOrbitsPos;
    orbits{i}{1} = orbitsPos;
    
    if showGraph
        [closedOrbitsNeg,orbitsNeg,hFigure(i,2)] = poincare_closed_orbit(domain,resolution,etaNeg,poincareSection,'odeSolverOptions',odeSolverOptions,'dThresh',dThresh,'nBisection',nBisection,'periodicBc',periodicBc,'showGraph',showGraph);
        hTitle = get(get(hFigure(i,2),'CurrentAxes'),'Title');
        titleString = get(hTitle,'String');
        newTitleString = sprintf('%s #%u \\eta_-',titleString,i);
        set(hTitle,'String',newTitleString)
    else
        [closedOrbitsNeg,orbitsNeg] = poincare_closed_orbit(domain,resolution,etaNeg,poincareSection,'odeSolverOptions',odeSolverOptions,'dThresh',dThresh,'nBisection',nBisection,'periodicBc',periodicBc,'showGraph',showGraph);
    end
    closedOrbits{i}{2} = closedOrbitsNeg;
    orbits{i}{2} = orbitsNeg;    
end

if showGraph
    varargout{1} = hFigure;
end
