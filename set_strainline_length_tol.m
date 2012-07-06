function strainline = set_strainline_length_tol(lengthTol,strainline)
% Set the strainline segment minimum length tolerance and delete fields
% that depend on it from the strainline structure.

validateattributes(lengthTol,{'double'},{'scalar','nonnegative'})

strainline.lengthTol = lengthTol;

fieldsToDelete = {'filteredSegmentIndex','hausdorffDistance',...
    'relativeStretching','segmentIndex'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
