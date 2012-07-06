function flow = set_flow_parameters(parameters,flow)

validateattributes(parameters,{'struct'},{})

flow.parameters = parameters;

fieldsToDelete = {'cgEigenvalue','cgEigenvector','cgStrain',...
    'finalPosition'};
    
for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
