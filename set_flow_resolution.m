% set_flow_resolution Set the flow resolution and delete fields that depend
% on it from the flow structure.
%
% SYNTAX
% flow = set_flow_resolution(resolution,flow)
%
% DESCRIPTION
% Set the flow resolution and delete fields that depend on it from the flow
% structure. Note that although the strainline and shearline structures
% have fields that depend on the flow structure fields, those fields in the
% strainline and shearline structures do not get delete by
% set_flow_resolution. The functions reset_strainline and reset_shearline
% can be used.
%
% INPUT ARGUMENTS
% resolution is a 1x2 array or a scalar. If it is a scalar, the number
% given is used for the x-resolution and the y-resolution is set
% automatically to have a square grid.
%
% EXAMPLE
% flow = set_flow_resolution([400 200],flow)
%
% SEE ALSO
% reset_strainline, reset_shearline

function flow = set_flow_resolution(resolution,flow)

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
