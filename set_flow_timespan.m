function output = set_flow_timespan(timespan,input)
% Set the flow timespan delete fields that depend on it from the flow,
% strainline and shearline structures.

validateattributes(timespan,{'double'},{'size',[1 2]})

if isfield(input,'flow')
    flow = input.flow;
    flow.timespan = timespan;
else
    error('Input has no flow field')
end

if isfield(flow,'finalPosition')
    flow = rmfield(flow,'finalPosition');
end

if isfield(flow,'cgStrain')
    flow = rmfield(flow,'cgStrain');
end

if isfield(flow,'cgEigenvector')
    flow = rmfield(flow,'cgEigenvector');
end

if isfield(flow,'cgEigenvalue')
    flow = rmfield(flow,'cgEigenvalue');
end

output.flow = flow;

if isfield(input,'strainline')
    % Use set_strainline_resolution to reset fields since it resets 
    % everything necessary.
    output.strainline = set_strainline_resolution(...
        input.strainline.resolution,input.strainline);
end

if isfield(input,'shearline')
    % Use set_shearline_resolution to reset fields since it resets 
    % everything necessary.
    output.shearline = set_shearline_resolution(...
        input.shearline.resolution,input.shearline);
end
