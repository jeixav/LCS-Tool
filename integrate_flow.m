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
else
    progressBar = [];
end

odeSolver = flow.odeSolver;
if ~isfield(flow,'odeSolverOptions')
    odeSolverOptions = [];
else
    odeSolverOptions = flow.odeSolverOptions;
end
odefun = flow.derivative;
timespan = flow.timespan;

if flow.coupledIntegration
    % FIXME coupledIntegration assumed equation of variation is used
    % FIXME Should not use for loop with coupled integration. This is a
    % quick fix.
    eovIc = [1 0 0 1];
    parfor iPosition = 1:nPosition
        flowSolution(iPosition) = feval(odeSolver,odefun,timespan,[initialPosition(iPosition,:) eovIc],odeSolverOptions);
        % Remove equation of variation terms
        flowSolution(iPosition).y = flowSolution(iPosition).y(1:2,:);
        % FIXME idata.f3d is not documented.
        flowSolution(iPosition).idata.f3d = flowSolution(iPosition).idata.f3d(1:2,:,:);
    end
else
    parfor iPosition = 1:nPosition
        flowSolution(iPosition) = feval(odeSolver,odefun,timespan,initialPosition(iPosition,:),odeSolverOptions);
        if verbose
            progressBar.increment(iPosition) %#ok<PFBNS>
        end
    end
end

if verbose
    try
        delete(progressBar)
    catch me %#ok<NASGU>
    end
end
