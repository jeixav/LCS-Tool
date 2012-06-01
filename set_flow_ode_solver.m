function output = set_flow_ode_solver(odeSolver,input)
% Set the flow integration method and delete fields that depend on it from
% the flow, strainline and shearline structures.

validateattributes(odeSolver,{'function_handle'},{})

if isfield(input,'flow')
    flow = input.flow;
else
    error('Input has no flow field')
end

switch func2str(odeSolver)
    case {'ode1','ode2','ode3','ode4','ode5','ode113','ode23','ode45',...
            'ode113v','ode23v','ode45v'}
        flow.odeSolver = odeSolver;
    otherwise
        error('ODE intgrator name not recognized.')
end
     
if isfield(flow,'odeSolverOptions')
    flow = rmfield(flow,'odeSolverOptions');
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
