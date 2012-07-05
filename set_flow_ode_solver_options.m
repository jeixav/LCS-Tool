function output = set_flow_ode_solver_options(odeSolverOptions,input)
% Set the flow timespan and delete fields that depend on it from the flow,
% strainline and shearline structures.

narginchk(1,2)
validateattributes(odeSolverOptions,{'struct'},{})

output = input;
clear('input')

if isfield(output,'flow')
    output.flow.odeSolverOptions = odeset(odeSolverOptions);
else
    error('Input has no flow field')
end

if isfield(output.flow,'finalPosition')
    output.flow = rmfield(output.flow,'finalPosition');
end

if isfield(output.flow,'cgStrain')
    output.flow = rmfield(output.flow,'cgStrain');
end

if isfield(output.flow,'cgEigenvector')
    output.flow = rmfield(output.flow,'cgEigenvector');
end

if isfield(output.flow,'cgEigenvalue')
    output.flow = rmfield(output.flow,'cgEigenvalue');
end

if isfield(output,'strainline')
    output.strainline = set_strainline_resolution(...
        output.strainline.resolution,output.strainline);
end

if isfield(output,'shearline')
    output.shearline = set_shearline_resolution(...
        output.shearline.resolution,output.shearline);
    if isfield(output.shearline,'etaPos')
        output.shearline = rmfield(output.shearline,'etaPos');
    end
    if isfield(output.shearline,'etaNeg')
        output.shearline = rmfield(output.shearline,'etaNeg');
    end

end
