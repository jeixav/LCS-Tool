function shearline = set_shearline_initial_position(initialPosition,...
    shearline)
% Set the shearline initial position(s) and delete fields that depend on 
% it from the shearline structure.

validateattributes(initialPosition,{'double'},{'size',[nan 2],'real',...
    'finite'})

shearline.initialPosition = initialPosition;

fieldsToDelete = {'averageGeodesicDeviationNeg',...
    'averageGeodesicDeviationPos','filteredIndexNeg','filteredIndexPos',...
    'geodesicDeviationNeg','geodesicDeviationPos','resolution',...
    'positionPos','positionNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
