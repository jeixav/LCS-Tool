%integrate_flow Integrate flow
%
% SYNTAX
% flowSolution = integrate_flow(flow,initialPosition,useEoV)
% flowSolution = integrate_flow(flow,initialPosition,useEoV,verbose)
%
% INPUT ARGUMENTS
% initialPosition: n-by-2 array
% useEoV: true or false
% verbose: true or false

function flowSolution = integrate_flow(flow,initialPosition,useEoV,verbose)

if ~isfield(flow,'odeSolver')
    flow.odeSolver = @ode45;
end

if nargin < 4
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
    initialPosition = transpose(initialPosition);
    initialPosition = initialPosition(:);
    flowSolution = feval(odeSolver,@(t,y)odefun(t,y,useEoV),timespan,initialPosition,odeSolverOptions);
else
    flowCgStrainMethodName = flow.cgStrainMethod.name;

    % FIXME If ODE solution structure is changed by MathWorks, this will
    % break
    flowSolution(nPosition) = struct('solver',[],'extdata',[],'x',[],'y',[],'stats',[],'idata',[]);

    switch flowCgStrainMethodName
        case 'finiteDifference'
            parfor iPosition = 1:nPosition
                flowSolution(iPosition) = feval(odeSolver,@(t,y)odefun(t,y,useEoV),timespan,initialPosition(iPosition,:),odeSolverOptions); %#ok<PFBNS>
                if verbose
                    progressBar.increment(iPosition) %#ok<PFBNS>
                end
            end
        case 'equationOfVariation'
            parfor iPosition = 1:nPosition
                if useEoV
                    dFlowMap0 = eye(2);
                    dFlowMap0 = reshape(dFlowMap0,4,1);
                    iInitialPosition = [initialPosition(iPosition,:),dFlowMap0];
                else
                    iInitialPosition = initialPosition(iPosition,:);
                end
                flowSolution(iPosition) = feval(odeSolver,@(t,y)odefun(t,y,useEoV),timespan,iInitialPosition,odeSolverOptions); %#ok<PFBNS>
                if verbose
                    progressBar.increment(iPosition) %#ok<PFBNS>
                end
                
            end
    end
end

if verbose
    try
        delete(progressBar)
    catch me %#ok<NASGU>
    end
end
