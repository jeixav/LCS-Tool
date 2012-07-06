function noStretchLine = set_no_stretch_line_resolution(resolution,...
    noStretchLine)
% Set the no stretch line initial positions resolution and delete fields
% that depend on it from the no stretch line structure.

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

noStretchLine.resolution = resolution;

fieldsToDelete = {'initialPosition','positionPos','positionNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(noStretchLine,fieldsToDelete{iField})
        noStretchLine = rmfield(noStretchLine,fieldsToDelete{iField});
    end
end
