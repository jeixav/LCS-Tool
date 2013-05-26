function flow = set_flow_domain(domain,flow)

validateattributes(domain,{'double'},{'size',[2 2]})

flow.domain = domain;

fieldsToDelete = {'cgEigenvalue','cgEigenvector','cgStrain','finalPosition','solution'};

for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
