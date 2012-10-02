function strainline = set_strainline_resolution(resolution,strainline)
% Set the strainline initial positions resolution and delete fields that
% depend on it from the strainline structure.

validateattributes(resolution,{'uint64'},{'size',[1 2],'positive'})

strainline.resolution = resolution;

fieldsToDelete = {'position','geodesicDeviation','segmentIndex',...
    'relativeStretching','hausdorffDistance','filteredSegmentIndex',...
    'initialPosition'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
