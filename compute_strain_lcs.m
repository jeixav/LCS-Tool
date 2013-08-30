function [flow,strainline] = compute_strain_lcs(flow,strainline,verbose)

narginchk(2,3)

verboseDefault = struct('progress',true,'stats',true);
if nargin < 3
    verbose = [];
end
verbose = set_default(verbose,verboseDefault);

% FIXME This if-statement is identical with one in compute_shear_lcs
if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    p = inputParser;
    p.KeepUnmatched = true;
    addParamValue(p,'cgStrainMethod',struct('name','finiteDifference'),@isstruct)
    addParamValue(p,'cgStrainCustomEigMethod',false,@islogical)
    addParamValue(p,'coupledIntegration',true,@islogical)
    parse(p,flow)
    
    [flow.cgEigenvalue,flow.cgEigenvector,flow.cgStrain] = eig_cgStrain(flow,p.Results.cgStrainMethod,p.Results.cgStrainCustomEigMethod,p.Results.coupledIntegration);
end

if ~isfield(strainline,'position')
    strainline = compute_strainline(flow,strainline,verbose);
end

if ~isfield(strainline,'geodesicDeviation')
    cgPosition = initialize_ic_grid(flow.resolution,flow.domain);
    strainline.geodesicDeviation = geodesic_deviation_strainline(strainline.position,cgPosition,flow.cgEigenvalue,flow.cgEigenvector,flow.resolution,verbose.progress);
    geodesic_deviation_stats(strainline.geodesicDeviation,true);
end

if ~isfield(strainline,'segmentIndex')
    strainline.segmentIndex = find_segments(strainline.position,strainline.geodesicDeviation,strainline.geodesicDeviationTol,strainline.lengthTol);
    nSegments = sum(cellfun(@(input)size(input,1),strainline.segmentIndex));
    disp(['Number of strainline segments: ',num2str(nSegments)])
end

if ~isfield(strainline,'relativeStretching')
    cgPosition = initialize_ic_grid(flow.resolution,flow.domain);
    strainline.relativeStretching = relative_stretching(strainline.position,strainline.segmentIndex,cgPosition,flow.cgEigenvalue(:,1),flow.resolution,verbose.progress);
end

if ~isfield(strainline,'filteredSegmentIndex')
    switch strainline.filteringMethod
        case 'superminimize'
            if ~isfield(verbose,'graphSuperminLine') || verbose.graphSuperminLine == false
                plotSuperminLine = false;
            else
                plotSuperminLine = verbose.graphSuperminLine;
            end
            matlabpoolClosed = false;
            if plotSuperminLine && matlabpool('size')
                warning([mfilename,':plotSuperminLine'],'plotSuperminLine does not work when matlabpool is in use. Temporarily closing matlabpool.')
                matlabpool('close')
                matlabpoolClosed = true;
            end
            strainline.filteredSegmentIndex = superminimize_grid(strainline.position,strainline.segmentIndex,strainline.relativeStretching,strainline.filteringParameters.distance,flow.domain,strainline.filteringParameters.resolution,plotSuperminLine,verbose.progress);
            if matlabpoolClosed
                matlabpool('open')
            end
        case 'hausdorff'
            strainline = hausdorff_filtering(strainline);
        case 'minimin'
            strainline = minimin_filtering(strainline);
        otherwise
            error('strainline.filteringMethod invalid')
    end
    nfilteredSegment = sum(cellfun(@sum,strainline.filteredSegmentIndex));
    disp(['Number of filtered segments: ' num2str(nfilteredSegment)])
end

function [segmentIndex,segmentLength] = find_segments(position,geodesicDeviation,geodesicDeviationTol,LengthTol)

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
    localLength = curve_length([position((candidateIndex(n,1):candidateIndex(n,2)),1),position((candidateIndex(n,1):candidateIndex(n,2)),2)]);
    if localLength >= LengthTol
        segmentIndex(n,:) = candidateIndex(n,:);
        segmentLength(n) = localLength;
    end
end
segmentIndex = segmentIndex(~isnan(segmentIndex(:,1)),:);
segmentLength = segmentLength(~isnan(segmentLength));

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

function length = curve_length(position)

x_d = diff(position(:,1));
y_d = diff(position(:,2));

length = sum(sqrt(x_d.^2 + y_d.^2));

function strainline = hausdorff_filtering(strainline,showPlot)

if nargin < 2
    showPlot = false;
end

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

% Hausdorff distance array. Stores the Hausdorff distance between all
% strainline segments.

if ~isfield(strainline,'hausdorffDistance')
    strainline.hausdorffDistance = nan(nSegments);
    % Set diagonal elements to 0
    strainline.hausdorffDistance(1:nSegments+1:end) = 0;

    strainlineSegmentPosition = strainline_segment_position(...
        strainline.position,strainline.segmentIndex);
        
    progressBar = ConsoleProgressBar;
    progressBar.setText('Hausdorff distance')
    progressBar.setTextPosition('left')
    progressBar.setElapsedTimeVisible(1);
    progressBar.setRemainedTimeVisible(1);
    progressBar.setLength(20)
    progressBar.setMaximum(nSegments);
    progressBar.start
    
    for iRefSegment = 1:nSegments
        iCompSegment = [1:(iRefSegment-1) (iRefSegment+1):nSegments];
        
        strainline.hausdorffDistance(iRefSegment,iCompSegment) = ...
            cellfun(@(input)hausdorff_dist(...
            strainlineSegmentPosition{iRefSegment},input),...
            strainlineSegmentPosition(iCompSegment));

        progressBar.setValue(progressBar.value + 1);
    end
    
    progressBar.setValue(progressBar.maximum)
    progressBar.stop
    fprintf('\n')
    
end

if showPlot
    hausdorffDistanceFigure = figure;

    hausdorffDistanceAxes = axes('nextplot','add',...
        'box','on',...
        'xgrid','on',...
        'ygrid','on',...
        'parent',hausdorffDistanceFigure);

    xlabel('Hausdorff distance','Parent',hausdorffDistanceAxes)
    ylabel('Relative stretching','Parent',hausdorffDistanceAxes)
end

for iStrainline = 1:nStrainlines
    nSegments = size(strainline.segmentIndex{iStrainline},1);
    for iSegment = 1:nSegments

        refIndex = (segmentIndex2(:,1) == iStrainline) & (segmentIndex2(:,2) == iSegment);
        lHausdorffDistance = strainline.hausdorffDistance(refIndex,:);
        minIndex = hausdorff_minimization([lHausdorffDistance; ...
            relativeStretchingLocal2].',strainline.filteringDistanceTol);

        if ~isempty(minIndex)
            strainline.filteredSegmentIndex{segmentIndex2(minIndex,1)}...
                (segmentIndex2(minIndex,2)) = true;
        end
            
        if showPlot
            % Kludge: Pass handle of main axes in showPlot
            strainlineAxes = showPlot;
            hWithinTol = plot_segment_within_tol(...
                strainline.filteringDistanceTol,lHausdorffDistance,...
                segmentIndex2,strainline.position,...
                strainline.segmentIndex,strainlineAxes);
            
            startIndex = strainline.segmentIndex{segmentIndex2(refIndex,1)}...
                (segmentIndex2(refIndex,2),1);
            stopIndex = strainline.segmentIndex{segmentIndex2(refIndex,1)}...
                (segmentIndex2(refIndex,2),2);
            refSegment = strainline.position{segmentIndex2(refIndex,1)}...
                (startIndex:stopIndex,:);
        
            hRefSegment = plot(strainlineAxes,...
                refSegment(:,1),refSegment(:,2),...
                'color','r','linewidth',1.5);

            p1 = plot(hausdorffDistanceAxes,lHausdorffDistance,...
                relativeStretchingLocal2,'LineStyle','none','Marker','o');
            if ~isempty(minIndex)
                p2 = plot(hausdorffDistanceAxes,...
                    lHausdorffDistance(minIndex),...
                    relativeStretchingLocal2(minIndex),'LineStyle','none',...
                    'Marker','o','MarkerFaceColor','r',...
                    'MarkerEdgeColor','k');
            end
            pause
            delete([p1; hWithinTol; hRefSegment])
            if ~isempty(minIndex)
                delete(p2)
            end
        end
    end
end

if showPlot
    delete(hausdorffDistanceFigure)
end

function minimumIndex = hausdorff_minimization(data,distanceTol)

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

function hSegment = plot_segment_within_tol(distanceTol,distance,...
    strainlineSegmentIndex,strainlinePosition,segmentIndex,strainlineAxes)

segmentWithinTol = strainlineSegmentIndex(distance < distanceTol,:);
nSegmentWithinTol = size(segmentWithinTol,1);

hSegment = nan(nSegmentWithinTol,1);

for iWithinTol = 1:nSegmentWithinTol
    
    iStrainline = segmentWithinTol(iWithinTol,1);
    iSegment = segmentWithinTol(iWithinTol,2);
    
    startIndex = segmentIndex{iStrainline}(iSegment,1);
    stopIndex = segmentIndex{iStrainline}(iSegment,2);
    
    position = strainlinePosition{iStrainline}(startIndex:stopIndex,:);

    hSegment(iWithinTol) = plot(strainlineAxes,position(:,1),position(:,2),'Color','k');
    
end

function [hd,D] = hausdorff_dist(P,Q)
% Calculates the Hausdorff Distance between P and Q
%
% hd = HausdorffDist(P,Q)
% [hd,D] = HausdorffDist(P,Q)
%
% Calculates the Hausdorff Distance, hd, between two sets of points, P and Q
% (which could be two trajectories) in two dimensions. Sets P and Q must,
% therefore, be matricies with an equal number of columns (dimensions),
% though not necessarily an equal number of rows (observations).
%
% The Directional Hausdorff Distance (dhd) is defined as:
% dhd(P,Q) = max p c P [ min q c Q [ ||p-q|| ] ].
% Intuitively dhd finds the point p from the set P that is farthest from any
% point in Q and measures the distance from p to its nearest neighbor in Q.
% 
% The Hausdorff Distance is defined as max{dhd(P,Q),dhd(Q,P)}
%
% D is the matrix of distances where D(n,m) is the distance of the nth
% point in P from the mth point in Q.
%
%
% 
% %%% ZCD Oct 2009 %%%
% Edits ZCD: Added the matrix of distances as an output. Fixed bug that
%   would cause an error if one of the sets was a single point. Removed
%   excess calls to "size" and "length". - May 2010
% Edits ZCD: Allowed for comparisons of N-dimensions. - June 2010
%

sP = size(P); sQ = size(Q);

if ~(sP(2)==sQ(2))
    error('Inputs P and Q must have the same number of columns')
end

% obtain all possible point comparisons
iP = repmat(1:sP(1),[1,sQ(1)])';
iQ = repmat(1:sQ(1),[sP(1),1]);
combos = [iP,iQ(:)];

% get distances for each point combination
cP=P(combos(:,1),:); cQ=Q(combos(:,2),:);
dists = sqrt(sum((cP - cQ).^2,2));

% Now create a matrix of distances where D(n,m) is the distance of the nth
% point in P from the mth point in Q. The maximum distance from any point
% in Q from P will be max(D,[],1) and the maximum distance from any point
% in P from Q will be max(D,[],1);
D = reshape(dists,sP(1),[]);

% Obtain the value of the point, p, in P with the largest minimum distance
% to any point in Q.
vp = max(min(D,[],2));
% Obtain the value of the point, q, in Q with the largets minimum distance
% to any point in P.
% vq = max(min(D,[],1));

% hd = max(vp,vq);
hd = vp;

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
