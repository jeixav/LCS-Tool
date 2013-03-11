% set_shearline_initial_position Set shearline initial positions
%
% SYNTAX
% shearline = set_shearline_initial_position(initialPosition,shearline)
% 
% DESCRIPTION
% Set shearline initial positions and delete fields that depend on them
% from the shearline structure.
%
% INPUT ARGUMENTS
% initialPosition: n-by-2 array
% shearline: LCS Toolbox shearline structure

function shearline = set_shearline_initial_position(initialPosition,shearline)

validateattributes(initialPosition,{'double'},{'size',[nan 2],'real','finite'})

shearline.initialPosition = initialPosition;

fieldsToDelete = {'averageGeodesicDeviationNeg','averageGeodesicDeviationPos','filteredIndexNeg','filteredIndexPos','geodesicDeviationNeg','geodesicDeviationPos','resolution','positionPos','positionNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
