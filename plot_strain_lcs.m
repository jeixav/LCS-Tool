function plot_strain_lcs(axes,flow,strainline,showPlot)

if isfield(flow,'cgEigenvector')
    cgPosition = initialize_ic_grid(flow.resolution,flow.domain);
    hQuiverPos = quiver(axes,cgPosition(:,1),cgPosition(:,2),flow.cgEigenvector(:,1),flow.cgEigenvector(:,2));
    set(hQuiverPos,'AutoScaleFactor',.25)
    set(hQuiverPos,'tag','quiver')
    set(hQuiverPos,'color','b')
    % Add negative eigenvectors to get double-headed arrows
    hQuiverNeg = quiver(axes,cgPosition(:,1),cgPosition(:,2),-flow.cgEigenvector(:,1),-flow.cgEigenvector(:,2));
    set(hQuiverNeg,'AutoScaleFactor',.25)
    set(hQuiverNeg,'tag','quiver')
    set(hQuiverNeg,'color','b')
    hQuiver = [hQuiverPos,hQuiverNeg];
    if ~isfield(showPlot,'quiver') || showPlot.quiver == false;
        set(hQuiver,'visible','off')
    end
end

if isfield(strainline,'position')
    hStrainline = cellfun(@(position)plot(axes,position(:,1),position(:,2)),strainline.position);
    set(hStrainline,'tag','strainline')
    set(hStrainline,'color',.8*[1 1 1])
    if ~isfield(showPlot,'strainline') || showPlot.strainline == false
        set(hStrainline,'visible','off')
    end
end

% FIXME May need to account for case when strainline.resolution is defined
% but strainline.initialPosition is not.
if isfield(strainline,'initialPosition')
    plot(axes,...
        strainline.initialPosition(:,1),strainline.initialPosition(:,2),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','o',...
        'LineStyle','none','Tag','strainlineInitialCondition')
    if ~isfield(showPlot,'strainlineInitialCondition') || ...
            showPlot.strainlineInitialCondition == false
        set(findobj(axes,'tag','strainlineInitialCondition'),'visible',...
            'off')
    end
end

if isfield(strainline,'geodesicDeviation')
    geodesicDeviationIndex = cellfun(...
        @(input)input<strainline.geodesicDeviationTol,...
        strainline.geodesicDeviation,'UniformOutput',false);
    cellfun(@(position,index)plot_geodesic_deviation_points(position,...
        index,axes),strainline.position,geodesicDeviationIndex)
    if ~isfield(showPlot,'geodesicDeviation') || ... 
            showPlot.geodesicDeviation == false
        set(findobj(axes,'tag','geodesicDeviationPoint'),'visible','off')
    end
end

if all(isfield(strainline,{'position','segmentIndex'}))
    plot_strainline_segment(axes,strainline.position,strainline.segmentIndex)
    if ~isfield(showPlot,'strainlineSegment') ||  showPlot.strainlineSegment == false
        set(findobj(axes,'tag','strainlineSegment'),'visible','off')
    end
end

if all(isfield(strainline,{'position','segmentIndex','filteredSegmentIndex'}))
    hFilteredStrainline = plot_filtered_strainline(axes,strainline.position,strainline.segmentIndex,strainline.filteredSegmentIndex);
    % Color repelling strainlines red
    if diff(flow.timespan) > 0
        set(hFilteredStrainline,'color','r')
    end
    if ~isfield(showPlot,'strainlineFiltered') || showPlot.strainlineFiltered == false
        set(hFilteredStrainline,'visible','off')
    end
end

if any(strcmp(strainline.filteringMethod,{'superminimize','superminimize2'}))
    plot_superminimization_lines(axes,flow.domain,...
        strainline.filteringParameters.resolution)
    if ~isfield(showPlot,'superminimizationLine') || ...
            showPlot.superminimizationLine == false
        set(findobj(axes,'tag','superminimizationLine'),'visible',...
            'off')
    end
end

function plot_geodesic_deviation_points(position,index,axes)
    
plot(axes,position(index,1),position(index,2),...
    'MarkerEdgeColor','r',...
    'Marker','o',...
    'LineStyle','none',...
    'Tag','geodesicDeviationPoint')

function plot_superminimization_lines(axes,flowDomain,superminResolution)

superminGrid = initialize_ic_grid(superminResolution,flowDomain);
superminGridX = reshape(superminGrid(:,1),fliplr(superminResolution));

for idx = 1:superminResolution(1)
    plot(axes,[superminGridX(1,idx) superminGridX(1,idx)],...
        [flowDomain(2,1) flowDomain(2,2)],...
        'linestyle','--',...
        'tag','superminimizationLine')
end

superminGridY = reshape(superminGrid(:,2),fliplr(superminResolution));
for idx = 1:superminResolution(2)
    plot(axes,[flowDomain(1,1) flowDomain(1,2)],...
        [superminGridY(idx,1) superminGridY(idx,1)],...
        'linestyle','--',...
        'tag','superminimizationLine')
end
