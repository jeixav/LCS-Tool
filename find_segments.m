function [segmentIndex segmentLength] = find_segments(position,...
    geodesicDeviation,geodesicDeviationTol,LengthTol)

pointIndex = cellfun(@(input)input<geodesicDeviationTol,...
    geodesicDeviation,'UniformOutput',false);

[segmentIndex,segmentLength] = cellfun(@(pointIndexArray,positionArray)...
    find_segments_array(pointIndexArray,positionArray,...
    LengthTol),pointIndex,position,'UniformOutput',false);

function [segmentIndex,segmentLength] = find_segments_array(pointIndex,...
    position,LengthTol)

candidateIndex = find_ones(transpose(pointIndex));

segmentIndex = nan(size(candidateIndex));
segmentLength = nan(size(candidateIndex,1),1);
for n = 1:size(candidateIndex,1)
    localLength = curve_length([...
        position((candidateIndex(n,1):candidateIndex(n,2)),1),...
        position((candidateIndex(n,1):candidateIndex(n,2)),2)]);
    if localLength >= LengthTol
        segmentIndex(n,:) = candidateIndex(n,:);
        segmentLength(n) = localLength;
    end
end
segmentIndex = segmentIndex(~isnan(segmentIndex(:,1)),:);
segmentLength = segmentLength(~isnan(segmentLength));
