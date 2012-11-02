% Set the strainline initial positions resolution and delete fields that
% depend on it from the strainline structure.

function strainline = set_strainline_resolution(resolution,strainline)

validateattributes(resolution,{'numeric'},{'size',[1 2],'positive'})

if ~isa(resolution,'uint64')
    resolution = uint64(resolution);
end

strainline.resolution = resolution;

fieldsToDelete = {'position','geodesicDeviation','segmentIndex',...
    'relativeStretching','hausdorffDistance','filteredSegmentIndex',...
    'initialPosition'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
