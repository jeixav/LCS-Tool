function closedLambdaLine = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,varargin)

p = inputParser;

% FIXME Add validationFcns common with eig_cgStrain, etc
addRequired(p,'domain')
addRequired(p,'resolution')
addRequired(p,'cgEigenvector')
addRequired(p,'cgEigenvalue')
addRequired(p,'lambda')
addRequired(p,'poincareSection')
addParameter(p,'forceEtaComplexNaN',false,@(i)validateattributes(i,{'logical'},{'scalar'}))
addParameter(p,'lambdaLineOdeSolverOptions',odeset)
addParameter(p,'showPoincareGraph',false,@(i)validateattributes(i,{'logical'},{'scalar'}))

parse(p,domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,varargin{:})

cgEigenvector = p.Results.cgEigenvector;
cgEigenvalue = p.Results.cgEigenvalue;
lambda = p.Results.lambda;
forceEtaComplexNaN = p.Results.forceEtaComplexNaN;
lambdaLineOdeSolverOptions = p.Results.lambdaLineOdeSolverOptions;
showPoincareGraph = p.Results.showPoincareGraph;

nPoincareSection = numel(poincareSection);
closedLambdaLineArea = zeros(1,nPoincareSection);
orbitArea = nan(1,2);
closedLambdaLine = cell(1,nPoincareSection);

for iLambda = lambda
    [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,iLambda,'forceComplexNaN',forceEtaComplexNaN);
    if showPoincareGraph
        [closedLambdaLineCandidate,~,hPoincareMap] = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',showPoincareGraph);
        for i = 1:nPoincareSection
            for j = 1:2 % etaPos,etaNeg
                hTitle = get(get(hPoincareMap(i,j),'CurrentAxes'),'Title');
                titleString = get(hTitle,'String');
                newTitleString = [titleString,' \lambda = ',num2str(iLambda)];
                set(hTitle,'String',newTitleString)
            end
        end
    else
        closedLambdaLineCandidate = poincare_closed_orbit_multi(domain,resolution,shearline,poincareSection,'odeSolverOptions',lambdaLineOdeSolverOptions,'showGraph',showPoincareGraph);
    end
        
    % keep outermost closed orbit
    for i = 1:nPoincareSection
        for j = 1:2 % etaPos,etaNeg
            orbitArea(j) = polyarea(closedLambdaLineCandidate{i}{j}{end}(:,1),closedLambdaLineCandidate{i}{j}{end}(:,2));
        end
        if max(orbitArea) > closedLambdaLineArea(i)
            closedLambdaLineArea(i) = max(orbitArea);
            closedLambdaLine{i}{1}{1} = closedLambdaLineCandidate{i}{1}{end};
            closedLambdaLine{i}{2}{1} = closedLambdaLineCandidate{i}{2}{end};
        end
    end
end
