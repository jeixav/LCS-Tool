function shearline = set_shearline_final_time(finalTime,shearline)
% Set the shearline integration final time and delete fields that depend on
% it from the shearline structure.

validateattributes(finalTime,{'double'},{'scalar','positive'})

shearline.finalTime = finalTime;

fieldsToDelete = {'averageGeodesicDeviationNeg',...
    'averageGeodesicDeviationPos','filteredIndexNeg','filteredIndexPos',...
    'geodesicDeviationNeg','geodesicDeviationPos','positionPos',...
    'positionNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
