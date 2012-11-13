% Set the flow integration method and delete fields that depend on it from
% the flow structure.

function flow = set_flow_ode_solver(odeSolver,flow)

validateattributes(odeSolver,{'function_handle'},{})

flow.odeSolver = odeSolver;

fieldsToDelete = {'cgEigenvalue','cgEigenvector','cgStrain',...
    'finalPosition'};

for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
