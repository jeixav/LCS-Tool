% set_no_stretch_line_final_time Set no stretch line final time
%
% DESCRIPTION
% noStretchLine = set_no_stretch_line_max_length(maxLength,noStretchLine)
%
% Set no stretch line maximum length. This length bounds integration times
% and is found heuristically.
%
% Fields in the shearline structure that depend on the maximum length
% are deleted.

function noStretchLine = set_no_stretch_line_max_length(maxLength,...
    noStretchLine)

validateattributes(maxLength,{'double'},{'scalar','positive'})

noStretchLine.maxLength = maxLength;

fieldsToDelete = {'positionNeg','positionPos'};

for iField = 1:length(fieldsToDelete)
    if isfield(noStretchLine,fieldsToDelete{iField})
        noStretchLine = rmfield(noStretchLine,fieldsToDelete{iField});
    end
end
