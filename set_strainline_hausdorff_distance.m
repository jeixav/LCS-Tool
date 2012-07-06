function strainline = set_strainline_hausdorff_distance(...
    hausdorffDistance,strainline)
% Set the Hausdorff distance tolerance and delete fields that depend on it
% from the strainline structure.

validateattributes(hausdorffDistance,{'double'},{'scalar','nonnegative'})

strainline.filteringDistanceTol = hausdorffDistance;

if isfield(strainline,'filteredSegmentIndex')
    strainline = rmfield(strainline,'filteredSegmentIndex');
end

