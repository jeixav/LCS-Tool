function strainline = set_strainline_length_tol(lengthTol,strainline)
% Set the strainline segment minimum length tolerance and delete fields
% that depend on it from the strainline structure.

validateattributes(lengthTol,{'double'},{'scalar','nonnegative'})

strainline.lengthTol = lengthTol;

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
