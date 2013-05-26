function flow = set_flow_derivative(derivative,varargin)

narginchk(1,2)

p = inputParser;
addRequired(p,'derivative',@(i)isa(i,'function_handle'));
addOptional(p,'flow',[],@isstruct);
parse(p,derivative,varargin{:})

flow = p.Results.flow;
flow.derivative = p.Results.derivative;

fieldsToDelete = {'cgEigenvalue','cgEigenvector','cgStrain','finalPosition','solution'};

for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
