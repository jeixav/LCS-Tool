% set_strainline_initial_position Set strainline initial positions
%
% SYNTAX
% strainline = set_strainline_initial_position(initialPosition,strainline)
% 
% DESCRIPTION
% Set strainline initial positions and delete fields that depend on them
% from the strainline structure.
%
% INPUT ARGUMENTS
% initialPosition: n-by-2 array
% strainline: LCS Toolbox strainline structure

function strainline = set_strainline_initial_position(initialPosition,strainline)

validateattributes(initialPosition,{'double'},{'size',[nan 2],'real','finite'})

strainline.initialPosition = initialPosition;

fieldsToDelete = {'filteredSegmentIndex','geodesicDeviation','position','relativeStretching','segmentIndex'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
