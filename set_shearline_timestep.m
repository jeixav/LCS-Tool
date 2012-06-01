function shearline = set_shearline_timestep(timestep,shearline)
% Set the shearline integration timestep and delete fields that depend on
% it from the shearline structure.

validateattributes(timestep,{'double'},{'scalar','positive'})

shearline.odeSolverOptions = odeset(shearline.odeSolverOptions,...
    'initialStep',timestep);

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

if isfield(shearline,'positionPos')
    shearline = rmfield(shearline,'positionPos');
end

if isfield(shearline,'positionNeg')
    shearline = rmfield(shearline,'positionNeg');
end
