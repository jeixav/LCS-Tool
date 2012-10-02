function flow = set_flow_timespan(timespan,flow)

validateattributes(timespan,{'double'},{'size',[1 2]})

flow.timespan = timespan;

fieldsToDelete = {'cgEigenvalue','cgEigenvector','cgStrain',...
    'finalPosition','solution'};

for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
