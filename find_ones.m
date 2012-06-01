function index = find_ones(vector)
%FIND_ONES Find consecutive ones in vector
%   Given a vector or ones and zeros, return indices to consecutive ones.
%   For example, given [0 1 1 0], return [2 3]. Given [0 1 1 0 1 1],
%   return {[2 3],[5 6]}.

% FIXME Try to rewrite this function using CIRCSHIFT

if ~isa(vector,'logical')
    error('find_ones:not_logical','Input vector not a logical array')
end

if size(vector,1) ~= 1
    error('find_ones:not_vector','Input array not a vector')
end

segments_found = 0;
next_position_to_check = 1;
% FIXME ceil(length(vector)/3) is probably more than the maximum possible
% number of segments
index = nan(ceil(length(vector)/3),2);
while next_position_to_check < length(vector)
    start_candidate = find(vector(next_position_to_check:end),1) + ...
        next_position_to_check - 1;
    if isempty(start_candidate)
        next_position_to_check = length(vector);
    else
        next_position_to_check = start_candidate+1;
    end
    if next_position_to_check <= length(vector)
        end_candidate = find(~vector(next_position_to_check:end),1)...
            +next_position_to_check-2;
        if isempty(end_candidate)
            end_candidate = length(vector);
        end
        if end_candidate ~= start_candidate
            index(segments_found+1,:) = [start_candidate end_candidate];
            segments_found = segments_found + 1;
            next_position_to_check = end_candidate + 2;
        else
            next_position_to_check = end_candidate + 1;
        end
    end
end

% http://www.mathworks.com/help/techdoc/data_analysis/f0-10104.html#f0-8511
index(any(isnan(index),2),:) = [];
