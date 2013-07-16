function [closedOrbits,orbits] = poincare_closed_orbit_multi(flow,shearline,PSList,odeSolverOptions,dThresh,showGraph)
% poincare_closed_orbit_multi(flow,shearline,PSList,odeSolverOptions) takes
% a list of specified poincare sections PSList and finds closed orbits
% around these poincare sections in the eta+ and eta- field given by shearline.
%
% INPUT
% PSList                                List of N poincare sections
% Format of PSList
% PSList{i}.endPosition = [p1x p1y; p2x p2y];
% PSList{i}.numPoints = nPoints;
% PSList{i}.integrationLength = [0 intLength];
% showGraph                             logical, set to 1 to show plots of poincare sections
%
% OUTPUT
% closedOrbits{}{}                      Positions of closed orbits
% Format of closeOrbits
% closedOrbits{i}{1}                    outermost closed orbit around poincare section i in eta+ field
% closedOrbits{i}{2}                    outermost closed orbit around poincare section i in eta- field
%
% orbits{}{}{}                          Positions of all orbits
%                                       Format: orbits{1}{2}{3}: 3rd {3}
%                                       orbit of 1st {1} poincare section
%                                       in eta- {2} field

narginchk(5,6)
if nargin == 5
    showGraph = true;
end

nPoincareSection = size(PSList,2);
closedOrbits = cell(1,nPoincareSection);
orbits = cell(1,nPoincareSection);

for i=1:nPoincareSection    
    fprintf('Searching closed orbits around poincare section %d ...\n', i);
    
    % define current poincare section
    poincareSection.endPosition = PSList{i}.endPosition;
    poincareSection.numPoints = PSList{i}.numPoints;
    poincareSection.integrationLength = PSList{i}.integrationLength;
    
    % etaPos
    etaField = shearline.etaPos;
    % find outermost orbit of each pointcare section
    [closedOrbitsPos, orbitsPos] = poincare_closed_orbit_mod(flow,...
        etaField,poincareSection,odeSolverOptions,5,dThresh,showGraph);
    closedOrbits{i}{1} = closedOrbitsPos;
    orbits{i}{1}       = orbitsPos;
    
    % etaNeg
    etaField = shearline.etaNeg;
    % find outermost orbit of each pointcare section
    [closedOrbitsNeg, orbitsNeg] = poincare_closed_orbit_mod(flow,...
        etaField,poincareSection,odeSolverOptions,5,dThresh,showGraph);
    closedOrbits{i}{2} = closedOrbitsNeg;
    orbits{i}{2}       = orbitsNeg;    
end

end
