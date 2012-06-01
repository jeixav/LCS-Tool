function output = load_input_file(filename)

[inputPath,name,ext] = fileparts(filename);

if strcmp(ext,'.m')
    previousPath = addpath(inputPath);
    output = feval(name);
    path(previousPath)
elseif strcmp(ext,'.mat')
    input = load(filename);
    names = fieldnames(input);
    output.flow = input.(names{1}).flow;
    output.shearline = input.(names{1}).shearline;
    output.strainline = input.(names{1}).strainline;
else
    error(['File extension not recognized (',ext,')'])
end

end