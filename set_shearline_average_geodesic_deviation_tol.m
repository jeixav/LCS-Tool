function shearline = set_shearline_average_geodesic_deviation_tol(...
geodesicDeviationTol,shearline)
% Set the strainline average geodesic deviation tolerance and delete fields
% that depend on it from the shearline structure.

validateattributes(geodesicDeviationTol,{'double'},{'size',[1 2],...
    'nonnegative'})

shearline.averageGeodesicDeviationPosTol = geodesicDeviationTol(1);
shearline.averageGeodesicDeviationNegTol = geodesicDeviationTol(2);

if isfield(shearline,'filteredIndexPos')
    shearline = rmfield(shearline,'filteredIndexPos');
end

if isfield(shearline,'filteredIndexNeg')
    shearline = rmfield(shearline,'filteredIndexNeg');
end

