% set_shearline_max_length Set shearline maximum length
%
% DESCRIPTION
% shearline = set_shearline_max_length(maxLength,shearline)
%
% Set shearline maximum length. This length bounds integration times and
% is found heuristically.
%
% Fields in the shearline structure that depend on the maximum length
% are deleted.

function shearline = set_shearline_max_length(maxLength,shearline)

validateattributes(maxLength,{'double'},{'scalar','positive'})

shearline.maxLength = maxLength;

fieldsToDelete = {'averageGeodesicDeviationNeg',...
    'averageGeodesicDeviationPos','filteredIndexNeg','filteredIndexPos',...
    'geodesicDeviationNeg','geodesicDeviationPos','positionPos',...
    'positionNeg'};

for iField = 1:length(fieldsToDelete)
    if isfield(shearline,fieldsToDelete{iField})
        shearline = rmfield(shearline,fieldsToDelete{iField});
    end
end
