function strainlineSegmentPosition = strainline_segment_position(...
    strainlinePosition,strainlineSegmentIndex)
%STRAINLINE_SEGMENT_POSITION Produce array of strainline segment positions
%   Given strainline.position and strainline.segmentIndex, return
%   strainline.segmentPosition. A cell array of strainline segments
%   positions.

nSegments = sum(cellfun(@(input)size(input,1),strainlineSegmentIndex));
nStrainlines = size(strainlinePosition,2);
strainlineSegmentPosition = cell(1,nSegments);

iSegment = 1;
for iStrainline = 1:nStrainlines;
    for iSegmentLocal = 1:size(strainlineSegmentIndex{iStrainline},1)
        startIndex = strainlineSegmentIndex{iStrainline}(iSegmentLocal,1);
        stopIndex = strainlineSegmentIndex{iStrainline}(iSegmentLocal,2);
        strainlineSegmentPosition{iSegment} = ...
            strainlinePosition{iStrainline}(startIndex:stopIndex,:);
    iSegment = iSegment + 1;
    end
end
        
end
