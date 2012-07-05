function output = set_flow_parameters(parameters,input)
% Set the flow parameters and delete fields that depend on it.

validateattributes(parameters,{'struct'},{})

if isfield(input,'flow')
    flow = input.flow;
    flow.parameters = parameters;
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

if isfield(input,'noStretchLine')
    output.noStretchLine = set_no_stretch_line_resolution(...
        input.noStretchLine.resolution,input.noStretchLine);
    if isfield(input.noStretchLine,'chiPos')
        output.noStretchLine = rmfield(output.noStretchLine,'chiPos');
    end
    if isfield(input.noStretchLine,'chiNeg')
        output.noStretchLine = rmfield(output.noStretchLine,'chiNeg');
    end
end
