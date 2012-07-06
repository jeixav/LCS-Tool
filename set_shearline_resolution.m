function shearline = set_shearline_resolution(resolution,shearline)
% Set the shearline initial positions resolution and delete fields that
% depend on it from the shearline structure.

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

shearline.resolution = resolution;

fieldsToDelete = {'averageGeodesicDeviationNeg',...
    'averageGeodesicDeviationPos','filteredIndexNeg','filteredIndexPos',...
    'geodesicDeviationNeg','geodesicDeviationPos','initialPosition',...
    'positionPos','positionNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
