function noStretchLine = set_no_stretch_line_ode_solver_options(...
    odeSolverOptions,noStretchLine)
% Set the no stretch line ODE solver options delete fields that depend on
% it from the no stretch line structure.

noStretchLine.odeSolverOptions = odeset(odeSolverOptions);

if isfield(noStretchLine,'positionPos')
    noStretchLine = rmfield(noStretchLine,'positionPos');
end

if isfield(noStretchLine,'positionNeg')
    noStretchLine = rmfield(noStretchLine,'positionNeg');
end
