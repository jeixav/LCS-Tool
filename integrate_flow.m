function flowSolution = integrate_flow(flow,initialPosition,verbose)

if ~isfield(flow,'odeSolver')
    flow.odeSolver = @ode45;
end

if nargin < 3
    verbose = true;
end

nPosition = size(initialPosition,1);

if verbose
    if ~exist('ParforProgressStarter2','file')
        addpath('ParforProgress2')
    end
    progressBar = ParforProgressStarter2(mfilename,nPosition);
end

odeSolver = flow.odeSolver;
if ~isfield(flow,'odeSolverOptions')
    odeSolverOptions = [];
else
    odeSolverOptions = flow.odeSolverOptions;
end
odefun = flow.derivative;
timespan = flow.timespan;
parfor iPosition = 1:nPosition
    flowSolution(iPosition) = feval(odeSolver,odefun,...
        timespan,initialPosition(iPosition,:),odeSolverOptions);
    if verbose
        progressBar.increment(iPosition) %#ok<PFBNS>
    end
end

if verbose
    try
        delete(progressBar)
    catch me %#ok<NASGU>
    end
end
