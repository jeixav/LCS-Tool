function output = set_flow_ode_solver_options(odeSolverOptions,input)
% Set the flow timespan and delete fields that depend on it from the flow,
% strainline and shearline structures.

validateattributes(odeSolverOptions,{'struct'},{})

output = input;
clear('input')

if isfield(output,'flow')
    output.flow.odeSolverOptions = odeset(odeSolverOptions);
else
    error('Input has no flow field')
end

fieldsToDelete = {'finalPosition','cgStrain','cgEigenvector',...
    'cgEigenvalue','strainline'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        output.flow = rmfield(output.flow,fieldsToDelete{iField});
    end
end

if isfield(output,'shearline')

    output.shearline = set_shearline_resolution(...
        output.shearline.resolution,output.shearline);
    fieldsToDelete = {'etaNeg','etaPos'};

    for iField = 1:length(fieldsToDelete)
        if isfield(output.shearline,fieldsToDelete{iField})
            output.shearline = rmfield(output.shearline,...
                fieldsToDelete{iField});
        end
    end
    
end

if isfield(output,'noStretchLine')

    output.noStretchLine = set_no_stretch_line_resolution(...
        output.noStretchLine.resolution,output.noStretchLine);
    fieldsToDelete = {'chiNeg','chiPos'};

    for iField = 1:length(fieldsToDelete)
        if isfield(output.noStretchLine,fieldsToDelete{iField})
            output.noStretchLine = rmfield(output.noStretchLine,...
                fieldsToDelete{iField});
        end
    end
    
end
