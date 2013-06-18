function [nextLocalMax,nextPosition] = find_next_local_max(array,flagArray,distance)

%% Find all local maxima with a distance threshold
[localMax,position] = local_max2D_gridded(array,distance);
[localMax,sortIndex] = sort(localMax,'descend');
position = position(sortIndex,:);

%% Check if next local maximum overlaps with flagArray
idx = 1;
while idx <= numel(localMax)
    if ~flagArray(position(idx,1),position(idx,2))
        nextLocalMax = localMax(idx);
        nextPosition = position(idx,:);
        return
    end
    idx = idx + 1;
end

% No local maxima left; return empty arrays
nextLocalMax = [];
nextPosition = [];
