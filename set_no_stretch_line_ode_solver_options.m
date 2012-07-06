function noStretchLine = set_no_stretch_line_ode_solver_options(...
    odeSolverOptions,noStretchLine)
% Set the no stretch line ODE solver options delete fields that depend on
% it from the no stretch line structure.

validateattributes(odeSolverOptions,{'struct'},{})

noStretchLine.odeSolverOptions = odeset(odeSolverOptions);

fieldsToDelete = {'positionNeg','positionPos'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
