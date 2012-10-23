function flowSolution = integrate_flow2(flow,initialPosition,verbose)

if ~isfield(flow,'odeSolver')
    flow.odeSolver = @ode45;
end

if nargin < 3
    verbose = true;
end

method = 'parfor';

switch method
    case 'arrayfun'
        initialPosition = to_coupled(initialPosition);
        
        flowSolution = ode_ic_array(flow.odeSolver,flow.derivative,flow.timespan,...
            initialPosition);
    case 'parfor'

        nPosition = prod(double(flow.resolution));

        if verbose
            progressBar = ParforProgressStarter2(mfilename,nPosition);
        end
        
        odeSolver = flow.odeSolver;
        odefun = flow.derivative;
        timespan = flow.timespan;
        parfor iPosition = 1:nPosition
            flowSolution(iPosition) = feval(odeSolver,odefun,...
                timespan,initialPosition(iPosition,:));
            if verbose
                progressBar.increment(iPosition) %#ok<PFBNS>
            end
        end
        
        if verbose
            try
                delete(progressBar);
            catch me %#ok<NASGU>
            end
        end

    otherwise
        error('Invalid method selected')
end

function sol = ode_ic_array(odeSolver,odefun,tspan,y0Array,options)

if nargin < 5
  options = [];
end

nTrajectories = size(y0Array,1)/2;

sol = arrayfun(@(idx)...
    feval(odeSolver,odefun,tspan,y0Array([idx,nTrajectories+idx]),...
    options),1:nTrajectories);
