function flow = set_flow_resolution(resolution,flow)
% Set the flow resolution and delete fields that depend on it from the
% flow, strainline and shearline structures.
%
% Example: doubleGyre = set_flow_resolution(uint64([400 200]),doubleGyre)

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

flow.resolution = resolution;

fieldsToDelete = {'cgEigenvalue','cgEigenvector','cgStrain',...
    'finalPosition'};

for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
