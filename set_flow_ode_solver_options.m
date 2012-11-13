% Set the flow ODE solver options structure and delete fields that depend
% on it from the flow structures.

function flow = set_flow_ode_solver_options(odeSolverOptions,flow)

validateattributes(odeSolverOptions,{'struct'},{})

flow.odeSolverOptions = odeSolverOptions;

fieldsToDelete = {'finalPosition','cgStrain','cgEigenvector',...
    'cgEigenvalue'};

for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
