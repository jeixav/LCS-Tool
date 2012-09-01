function strainline = set_strainline_final_time(finalTime,strainline)
% Set the strainline initial positions resolution and delete fields that
% depend on it from the strainline structure.

validateattributes(finalTime,{'double'},{'scalar','positive'})

strainline.finalTime = finalTime;

fieldsToDelete = {'position','geodesicDeviation','segmentIndex',...
    'relativeStretching','hausdorffDistance','filteredSegmentIndex'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
