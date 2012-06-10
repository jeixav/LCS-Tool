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

end
