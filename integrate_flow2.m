function flowSolution = integrate_flow2(flow,initialPosition)

initialPosition = to_coupled(initialPosition);

if ~isfield(flow,'odeSolver')
    flow.odeSolver = @ode45;
end

flowSolution = ode_ic_array(flow.odeSolver,flow.derivative,flow.timespan,...
    initialPosition);

function sol = ode_ic_array(odeSolver,odefun,tspan,y0Array,options)

if nargin < 5
  options = [];
end

nTrajectories = size(y0Array,1)/2;

sol = arrayfun(@(idx)...
    feval(odeSolver,odefun,tspan,y0Array([idx,nTrajectories+idx]),...
    options),1:nTrajectories);
