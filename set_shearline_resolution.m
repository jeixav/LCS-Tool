function shearline = set_shearline_resolution(resolution,shearline)
% Set the strainline initial positions resolution and delete fields that
% depend on it from the strainline structure.

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

shearline.resolution = resolution;

if isfield(shearline,'averageGeodesicDeviationNeg')
    shearline = rmfield(shearline,'averageGeodesicDeviationNeg');
end

if isfield(shearline,'averageGeodesicDeviationPos')
    shearline = rmfield(shearline,'averageGeodesicDeviationPos');
end

if isfield(shearline,'filteredIndexNeg')
    shearline = rmfield(shearline,'filteredIndexNeg');
end

if isfield(shearline,'filteredIndexPos')
    shearline = rmfield(shearline,'filteredIndexPos');
end

if isfield(shearline,'geodesicDeviationNeg')
    shearline = rmfield(shearline,'geodesicDeviationNeg');
end

if isfield(shearline,'geodesicDeviationPos')
    shearline = rmfield(shearline,'geodesicDeviationPos');
end

if isfield(shearline,'initialPosition')
    shearline = rmfield(shearline,'initialPosition');
end

if isfield(shearline,'positionPos')
    shearline = rmfield(shearline,'positionPos');
end

if isfield(shearline,'positionNeg')
    shearline = rmfield(shearline,'positionNeg');
end
