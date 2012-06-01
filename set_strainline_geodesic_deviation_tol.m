function strainline = set_strainline_geodesic_deviation_tol(...
    geodesicDeviationTol,strainline)
% Set the strainline segment maximum geodesic deviation tolerance and
% delete fields that depend on it from the strainline structure.

validateattributes(geodesicDeviationTol,{'double'},{'scalar','nonnegative'})

strainline.geodesicDeviationTol = geodesicDeviationTol;

if isfield(strainline,'segmentIndex')
    strainline = rmfield(strainline,'segmentIndex');
end

if isfield(strainline,'relativeStretching')
    strainline = rmfield(strainline,'relativeStretching');
end

if isfield(strainline,'hausdorffDistance')
    strainline = rmfield(strainline,'hausdorffDistance');
end

if isfield(strainline,'filteredSegmentIndex')
    strainline = rmfield(strainline,'filteredSegmentIndex');
end

end
