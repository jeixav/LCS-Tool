% set_strainline_max_length Set strainline maximum length
%
% DESCRIPTION
% strainline = set_strainline_max_length(maxLength,strainline)
%
% Set strainline maximum length. Generally strainline integration will
% stop at domain boundaries using ODE event detection. Nonetheless a 
% maximum length needs to be set to bound integration times. This maximum 
% length is found heuristically.
%
% Fields in the strainline structure that depend on the maximum length
% are deleted.

function strainline = set_strainline_max_length(maxLength,strainline)
 
validateattributes(maxLength,{'double'},{'scalar','positive'})

strainline.maxLength = maxLength;

fieldsToDelete = {'position','geodesicDeviation','segmentIndex',...
    'relativeStretching','hausdorffDistance','filteredSegmentIndex'};

for iField = 1:length(fieldsToDelete)
    if isfield(strainline,fieldsToDelete{iField})
        strainline = rmfield(strainline,fieldsToDelete{iField});
    end
end
