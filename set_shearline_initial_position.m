function shearline = set_shearline_initial_position(initialPosition,...
    shearline)
% Set the shearline initial position(s) and delete fields that depend on 
% it from the shearline structure.

validateattributes(initialPosition,{'double'},{'size',[nan 2],'real','finite'})

shearline.initialPosition = initialPosition;

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

if isfield(shearline,'resolution')
    shearline = rmfield(shearline,'resolution');
end

if isfield(shearline,'positionPos')
    shearline = rmfield(shearline,'positionPos');
end

if isfield(shearline,'positionNeg')
    shearline = rmfield(shearline,'positionNeg');
end
