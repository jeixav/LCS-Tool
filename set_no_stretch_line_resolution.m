function noStretchLine = set_no_stretch_line_resolution(resolution,...
    noStretchLine)
% Set the no stretch line initial positions resolution and delete fields
% that depend on it from the no stretch line structure.

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

noStretchLine.resolution = resolution;

if isfield(noStretchLine,'initialPosition')
    noStretchLine = rmfield(noStretchLine,'initialPosition');
end

if isfield(noStretchLine,'positionPos')
    noStretchLine = rmfield(noStretchLine,'positionPos');
end

if isfield(noStretchLine,'positionNeg')
    noStretchLine = rmfield(noStretchLine,'positionNeg');
end
