function flow = set_flow_resolution(resolution,flow)
% Set the flow resolution and delete fields that depend on it from the
% flow, strainline and shearline structures.
%
% The resolution can be specified as a 1x2 matrix or a scalar. If it is a
% scalar, the number given is used for the x-resolution and the
% y-resolution is set automatically to have a square grid.
%
% Example: flow = set_flow_resolution([400 200],flow)

validateattributes(resolution,{'numeric'},{'row','positive'})

if isscalar(resolution)
    resX = resolution;
    resY = resX/diff(flow.domain(1,:))*diff(flow.domain(2,:));
    resolution = [resX resY];
end
    
if ~isa(resolution,'uint64')
    resolution = uint64(resolution);
end

flow.resolution = resolution;

fieldsToDelete = {'cgEigenvalue','cgEigenvector','cgStrain',...
    'finalPosition','solution'};

for iField = 1:length(fieldsToDelete)
    if isfield(flow,fieldsToDelete{iField})
        flow = rmfield(flow,fieldsToDelete{iField});
    end
end
