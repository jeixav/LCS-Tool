function strainline = minimin_filtering(strainline)

nStrainlines = size(strainline.position,2);

strainline.filteredSegmentIndex = arrayfun(...
    @(idx)false(size(strainline.segmentIndex{idx},1),1),...
    1:nStrainlines,...
    'UniformOutput',false);

nSegments = sum(cellfun(@(input)size(input,1),strainline.segmentIndex));

segmentIndex2 = nan(nSegments,2);
counter = 1;
for m = 1:size(strainline.segmentIndex,2)
    for n = 1:size(strainline.segmentIndex{m},1)
        segmentIndex2(counter,1) = m;
        segmentIndex2(counter,2) = n;
        counter = counter + 1;
    end
end

relativeStretchingLocal2 = nan(1,nSegments);
for m = 1:nSegments
    relativeStretchingLocal2(m) = ...
        strainline.relativeStretching{segmentIndex2(m,1)}(segmentIndex2(m,2));
end

if ~isfield(strainline,'miniminDistance')
    strainline.miniminDistance = nan(nSegments);
    % Set diagonal elements to 0
    strainline.miniminDistance(1:nSegments+1:end) = 0;

    strainlineSegmentPosition = strainline_segment_position(...
        strainline.position,strainline.segmentIndex);
        
    progressBar = ParforProgressStarter2(mfilename,nSegments);
    
    iCompSegment = cell(1,nSegments - 1);
    for iRefSegment = 1:nSegments - 1
        iCompSegment{iRefSegment} = (iRefSegment + 1):nSegments;
    end
    
    parfor iRefSegment = 1:nSegments - 1
%         iCompSegment = (iRefSegment + 1):nSegments;
        
%         strainline.miniminDistance(iRefSegment,iCompSegment) = ...
%             cellfun(@(input)minimin_dist(...
%             strainlineSegmentPosition{iRefSegment},input),...
%             strainlineSegmentPosition(iCompSegment));

        miniminDistance{iRefSegment} = ...
            cellfun(@(input)minimin_dist(...
            strainlineSegmentPosition{iRefSegment},input),...
            strainlineSegmentPosition(iCompSegment{iRefSegment})); %#ok<PFBNS>

        progressBar.increment(iRefSegment); %#ok<PFBNS>
    end
    
    try
        delete(progressBar);
    catch me %#ok<NASGU>
    end
    
    for iRefSegment = 1:nSegments - 1
        strainline.miniminDistance(iRefSegment,iCompSegment{iRefSegment}) ...
            = miniminDistance{iRefSegment};
    end
    strainline.miniminDistance = triu(strainline.miniminDistance) ...
        + tril(transpose(strainline.miniminDistance));
end

for iStrainline = 1:nStrainlines
    nSegments = size(strainline.segmentIndex{iStrainline},1);
    for iSegment = 1:nSegments

        refIndex = (segmentIndex2(:,1) == iStrainline) ...
            & (segmentIndex2(:,2) == iSegment);
        lMiniminDistance = strainline.miniminDistance(refIndex,:);
        minIndex = minimin_minimization([lMiniminDistance; ...
            relativeStretchingLocal2].',...
            strainline.filteringParameters.distance);

        if ~isempty(minIndex)
            strainline.filteredSegmentIndex{segmentIndex2(minIndex,1)}...
                (segmentIndex2(minIndex,2)) = true;
        end
    end
end

function minimumIndex = minimin_minimization(data,distanceTol)

distance = data(:,1);
metric = data(:,2);

withinDistanceTolIndex = distance <= distanceTol;

[~,minIndex] = min(metric(withinDistanceTolIndex));

tmp = find(withinDistanceTolIndex);
if distance(tmp(minIndex)) == 0
    minimumIndex = tmp(minIndex);
else
    minimumIndex = [];
end

% Calculates the minimin distance between P and Q
function d = minimin_dist(p,q)

sP = size(p);
sQ = size(q);

if ~(sP(2)==sQ(2))
    error('Inputs P and Q must have the same number of columns')
end

% obtain all possible point comparisons
iP = repmat(1:sP(1),[1,sQ(1)])';
iQ = repmat(1:sQ(1),[sP(1),1]);
combos = [iP,iQ(:)];

% get distances for each point combination
cP = p(combos(:,1),:);
cQ = q(combos(:,2),:);
dists = sqrt(sum((cP - cQ).^2,2));

% Now create a matrix of distances where D(n,m) is the distance of the nth
% point in P from the mth point in Q. The maximum distance from any point
% in Q from P will be max(D,[],1) and the maximum distance from any point
% in P from Q will be max(D,[],1);
d = reshape(dists,sP(1),[]);

d = min(d(:));
