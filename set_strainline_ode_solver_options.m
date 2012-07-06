function strainline = set_strainline_ode_solver_options(...
    odeSolverOptions,strainline)
% Set the strainline ODE solver options and delete fields that depend on it
% from the strainline structure.

validateattributes(odeSolverOptions,{'struct'},{})

strainline.odeSolverOptions = odeset(odeSolverOptions);

fieldsToDelete = {'filteredSegmentIndex','geodesicDeviation',...
    'hausdorffDistance','position','relativeStretching','segmentIndex'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
