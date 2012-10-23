function shearline = set_shearline_ode_solver_options(odeSolverOptions,...
    shearline)
% Set the shearline ODE solver options delete fields that depend on it
% from the shearline structure.

validateattributes(odeSolverOptions,{'struct'},{})

shearline.odeSolverOptions = odeset(odeSolverOptions);

fieldsToDelete = {'averageGeodesicDeviationNeg',...
    'averageGeodesicDeviationPos','filteredIndexNeg','filteredIndexPos',...
    'geodesicDeviationNeg','geodesicDeviationPos','positionPos',...
    'positionNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
