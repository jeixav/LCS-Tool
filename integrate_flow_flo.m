function finalPosition = integrate_flow_flo(flow,initialPosition,verbose)
% INPUT:        
% initialPosition       Nx2 array with initial positions
% OUTPUT:
% finalPosition         Nx2 array with final positions

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

% coupled integration or one intial conditions one after the other
if flow.coupledIntegration
    % reshape initial positions into one 2Nx1 vector
    initialPositionOneColumn(1:2:2*nPosition-1) = initialPosition(:,1);
    initialPositionOneColumn(2:2:2*nPosition) = initialPosition(:,2);
    flowSolution = feval(odeSolver,@(t,y)odefun(t,y,false),timespan,initialPositionOneColumn,odeSolverOptions);    
    finalPos = arrayfun(@(odeSolution)deval(odeSolution,...
        flow.timespan(2)),flowSolution,'uniformOutput',...
        false);
    finalPos = cell2mat(finalPos);
    clear finalPosition
    finalPosition(:,1) = finalPos(1:2:end-1);
    finalPosition(:,2) = finalPos(2:2:end);        
else
    parfor iPosition = 1:nPosition
        flowSolution(iPosition) = feval(odeSolver,odefun,timespan,initialPosition(iPosition,:),odeSolverOptions);
        if verbose
            progressBar.increment(iPosition) %#ok<PFBNS>
        end
    end    
    finalPosition = arrayfun(@(odeSolution)deval(odeSolution,...
        flow.timespan(2)),flowSolution,'uniformOutput',...
        false);
    finalPosition = cell2mat(finalPosition);
    finalPosition = transpose(finalPosition);
end

if verbose
    try
        delete(progressBar)
    catch me %#ok<NASGU>
    end
end
