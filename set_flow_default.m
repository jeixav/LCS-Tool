function flow = set_flow_default(flow)
% Add default values to flow structure.

if ~isfield(flow,'isCompressible')
    flow.isCompressible = true;
end

if ~isfield(flow,'odeSolver')
    flow.odeSolver = @ode45;
end

if ~isfield(flow,'odeSolverOptions')
    switch func2str(flow.odeSolver)
        case {'ode1','ode2','ode3','ode4','ode5'}
            defaultTimestep = sign(diff(flow.timespan))*.1;
            % Kludge: For fixed timestep integrators, store timestep in
            % odeset's InitialStep parameter.
            flow.odeSolverOptions = odeset('InitialStep',defaultTimestep);
        case {'ode113','ode23','ode45','ode113v','ode23v','ode45v'}
            flow.odeSolverOptions = [];
        otherwise
            error('ODE solver type not known.')
    end
end
