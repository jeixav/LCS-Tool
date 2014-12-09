% poincare_closed_orbit_range Find closed orbits over ranges of lambda
%
% SYNTAX 
% [closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection)
%
% OUTPUT ARGUMENTS
% closedLambdaLinePos: closed lambda line positions for etaPos, cell array
% of size(numel(lambda),numel(poincareSection))
%
% closedLambdaLineNeg: closed lambda line positions for etaPos, cell array
% of size(numel(lambda),numel(poincareSection))
%
% EXAMPLES
% If lambda = [.99,1,1.01] and poincareSection is a 1x2 struct, to plot the
% innermost closed orbit of lambda(2), poincareSection(1):
% plot(closedLambdaLinePos{2,1}{1}(:,1),closedLambdaLinePos{2,1}{1}(:,2))
%
% Area enclosed by the outermost orbit:
% polyarea(closedLambdaLinePos{2,1}{end}(:,1),closedLambdaLinePos{2,1}{end}(:,2))

function [closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,varargin)

p = inputParser;

% FIXME Add validationFcns common with eig_cgStrain, etc
addRequired(p,'domain')
addRequired(p,'resolution')
addRequired(p,'cgEigenvector')
addRequired(p,'cgEigenvalue')
addRequired(p,'lambda')
addRequired(p,'poincareSection')
addParameter(p,'forceEtaComplexNaN',false,@(i)validateattributes(i,{'logical'},{'scalar'}))
addParameter(p,'odeSolverOptions',odeset)
addParameter(p,'periodicBc',[false,false],@(input)validateattributes(input,{'logical'},{'size',[1,2]}));
addParameter(p,'showPoincareGraph',false,@(i)validateattributes(i,{'logical'},{'scalar'}))

parse(p,domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,varargin{:})

forceEtaComplexNaN = p.Results.forceEtaComplexNaN;
odeSolverOptions = p.Results.odeSolverOptions;
periodicBc = p.Results.periodicBc;
showPoincareGraph = p.Results.showPoincareGraph;

nLambda = numel(lambda);
nPoincareSection = numel(poincareSection);
closedLambdaLinePos = cell(nLambda,nPoincareSection);
closedLambdaLineNeg = cell(nLambda,nPoincareSection);

% eta = [etaPos,etaNeg]
nEta = 2;

for iLambda = 1:nLambda
    [etaPos,etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda(iLambda),'forceComplexNaN',forceEtaComplexNaN);
    if showPoincareGraph
        [closedLambdaLine,~,hPoincareMap] = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,poincareSection,'odeSolverOptions',odeSolverOptions,'periodicBc',periodicBc,'showGraph',showPoincareGraph);
        for iPoincareSection = 1:nPoincareSection
            for iEta = 1:nEta
                hTitle = get(get(hPoincareMap(iPoincareSection,iEta),'CurrentAxes'),'Title');
                titleString = get(hTitle,'String');
                newTitleString = [titleString,' \lambda = ',num2str(lambda(iLambda))];
                set(hTitle,'String',newTitleString)
            end
        end
    else
        closedLambdaLine = poincare_closed_orbit_multi(domain,resolution,etaPos,etaNeg,poincareSection,'odeSolverOptions',odeSolverOptions,'periodicBc',periodicBc);
    end
    for iPoincareSection = 1:nPoincareSection
        closedLambdaLinePos{iLambda,iPoincareSection} = closedLambdaLine{iPoincareSection}{1};
        closedLambdaLineNeg{iLambda,iPoincareSection} = closedLambdaLine{iPoincareSection}{2};
    end
end
