function strainline = set_strainline_resolution(resolution,strainline)
% Set the strainline initial positions resolution and delete fields that
% depend on it from the strainline structure.

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

strainline.resolution = resolution;

if isfield(strainline,'position')
    strainline = rmfield(strainline,'position');
end

if isfield(strainline,'geodesicDeviation')
    strainline = rmfield(strainline,'geodesicDeviation');
end

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
