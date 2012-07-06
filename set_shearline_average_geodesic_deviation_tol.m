function shearline = set_shearline_average_geodesic_deviation_tol(...
geodesicDeviationTol,shearline)
% Set the shearline average geodesic deviation tolerances and delete fields
% that depend on it from the shearline structure.

validateattributes(geodesicDeviationTol,{'double'},{'size',[1 2],...
    'nonnegative'})

shearline.averageGeodesicDeviationPosTol = geodesicDeviationTol(1);
shearline.averageGeodesicDeviationNegTol = geodesicDeviationTol(2);

fieldsToDelete = {'filteredIndexPos','filteredIndexNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
