function finalPosition = integrate_flow(flow,initialPosition)

% Derivative,timespan,initialPosition,odeSolver,odeSolverOptions)

switch size(initialPosition,2)
    case 2
        auxiliary_grid_flag = false;
    case 8
        auxiliary_grid_flag = true;
    otherwise
        error('Unable to determine if auxiliary grid values are present.')
end

initialPosition = to_coupled(initialPosition);

if ~isempty(odeget(flow.odeSolverOptions,'OutputFcn')) && ...
        strcmp(func2str(odeget(flow.odeSolverOptions,'OutputFcn')),...
        'ode_progress_bar')
    progressBar = true;
    cpb = ConsoleProgressBar;
    cpb.setText(mfilename)
    cpb.setTextPosition('left')
    cpb.setElapsedTimeVisible(1);
    cpb.setRemainedTimeVisible(1);
    cpb.setLength(20)
    cpb.start
else
    progressBar = false;
end

if isfield(flow,'parameters')
    flowDerivative = @(t,y)flow.derivative(t,y,flow.parameters);
elseif isfield(handles.flow,'data')
    flowDerivative = @(t,y)flow.derivative(t,y,flow.data);
else
    flowDerivative = @(t,y)flow.derivative(t,y);
end

switch func2str(flow.odeSolver)
    % case 'fixedTimestep'
    case {'ode1','ode2','ode3','ode4','ode5','ode113v','ode23v','ode45v'}
        fixedOdeSolverOptions = flow.odeSolverOptions;
        if progressBar
            cpb.setMaximum(abs(diff(flow.timespan)))
            fixedOdeSolverOptions.OutputFcn = @(t,y,flag)...
                coupled_progress_bar(t,y,flag,cpb,flow.timespan(1));
        end
        if any(strcmp(func2str(flow.odeSolver),{'ode1','ode2','ode3','ode4','ode5'}))
            timestep = odeget(flow.odeSolverOptions,'initialStep');
            timespan = flow.timespan(1):timestep:flow.timespan(end);
            % flowDerivativeVector = @(t,y)cell2mat(arrayfun(@(m)flowDerivative(t,y(m,:)),1:size(y,1),'uniformOutput',false));
            finalPosition = feval(flow.odeSolver,flowDerivative,timespan,...
                initialPosition,fixedOdeSolverOptions);
        end
        if any(strcmp(func2str(flow.odeSolver),{'ode113v','ode23v','ode45v'}))
            % Remove last character of flow.odeSolver to obtain name of 
            % MATLAB built-in ODE solver.
            odeSolver = func2str(flow.odeSolver);
            odeSolver = odeSolver(1:end-1);
            odeSolver = str2func(odeSolver);
            [~,finalPosition] = feval(odeSolver,flowDerivative,flow.timespan,...
                initialPosition,fixedOdeSolverOptions);
        end
        finalPosition = transpose(finalPosition(end,:));
    % case 'variableTimestep'
    case {'ode113','ode23','ode45'}
        variableOdeSolverOptions = flow.odeSolverOptions;
        if progressBar
            cpb.setMaximum(size(initialPosition,1)/2);
            variableOdeSolverOptions.OutputFcn = @(t,y,flag)...
                array_progress_bar(t,y,flag,cpb);
        end
        solution = ode_ic_array(flow.odeSolver,flowDerivative,flow.timespan,...
            initialPosition,variableOdeSolverOptions);
        finalPositionCell = arrayfun(@(iSolution)deval(iSolution,...
            flow.timespan(end)),solution,'UniformOutput',false);
        nTrajectories = size(initialPosition,1)/2;
        finalPosition = reshape(transpose(cell2mat(finalPositionCell)),...
            nTrajectories*2,1);
    otherwise
        error('ODE intgrator type not known')
end

if progressBar
    cpb.setValue(cpb.maximum)
    cpb.stop
    fprintf('\n')
end

finalPosition = from_coupled(finalPosition,auxiliary_grid_flag);

function sol = ode_ic_array(odeSolver,odefun,tspan,y0Array,options)

if nargin < 5
  options = [];
end

nTrajectories = size(y0Array,1)/2;

sol = arrayfun(@(idx)...
    feval(odeSolver,odefun,tspan,y0Array([idx,nTrajectories+idx]),...
    options),1:nTrajectories);

function status = coupled_progress_bar(t,~,flag,progressBar,initialTime)

if nargin < 3 || isempty(flag)
    progressBar.setValue(abs(t(end)-initialTime))
end

status = 0;

function status = array_progress_bar(~,~,flag,progressBar)

if ~isempty(flag) && strcmp(flag,'done')
    progressBar.setValue(progressBar.value + 1);
end

status = 0;
