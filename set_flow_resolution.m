function output = set_flow_resolution(resolution,input)
% Set the flow resolution and delete fields that depend on it from the
% flow, strainline and shearline structures.
%
% Example: doubleGyre = set_flow_resolution(uint64([400 200]),doubleGyre)

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

if isfield(input,'flow')
    flow = input.flow;
    flow.resolution = resolution;
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
    output.shearline = set_shearline_resolution(...
        input.shearline.resolution,input.shearline);
    if isfield(input.shearline,'etaPos')
        output.shearline = rmfield(output.shearline,'etaPos');
    end
    if isfield(input.shearline,'etaNeg')
        output.shearline = rmfield(output.shearline,'etaNeg');
    end
end
