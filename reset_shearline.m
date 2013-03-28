% Remove all calculated fields of a shearline structure.

function shearline = reset_shearline(shearline)

fieldsToDelete = {'etaPos','etaNeg','positionPos','positionNeg','geodesicDeviationPos','geodesicDeviationNeg','averageGeodesicDeviationPos','averageGeodesicDeviationNeg','filteredIndexPos','filteredIndexNeg','geodesicDeviationPosMeshgrid','geodesicDeviationNegMeshgrid'};
    
for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
