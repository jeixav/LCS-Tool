function strainline = set_strainline_geodesic_deviation_tol(...
    geodesicDeviationTol,strainline)
% Set the strainline segment maximum geodesic deviation tolerance and
% delete fields that depend on it from the strainline structure.

validateattributes(geodesicDeviationTol,{'double'},{'scalar','nonnegative'})

strainline.geodesicDeviationTol = geodesicDeviationTol;

fieldsToDelete = {'filteredSegmentIndex','hausdorffDistance',...
    'relativeStretching','segmentIndex'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end

