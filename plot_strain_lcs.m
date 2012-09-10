function plot_strain_lcs(axes,flow,strainline,showPlot)

if isfield(strainline,'cgEigenvector')
    cgPosition = initialize_ic_grid(flow.resolution,flow.domain);
    quiver(axes,cgPosition(:,1),cgPosition(:,2),...
        flow.cgEigenvector(:,1),flow.cgEigenvector(:,2),'tag','quiver')
    if ~isfield(showPlot,'quiver') || showPlot.quiver == false;
        set(findobj(axes,'tag','quiver'),'visible','off')
    end
end

if isfield(strainline,'position')
    plot_strainline(axes,strainline)
    if ~isfield(showPlot,'strainline') || showPlot.strainline == false
        set(findobj(axes,'tag','strainline'),'visible','off')
    end
end

strainlineIc = initialize_ic_grid(strainline.resolution,flow.domain);
plot(axes,strainlineIc(:,1),strainlineIc(:,2),...
    'MarkerFaceColor','k',...
    'MarkerEdgeColor','k',...
    'Marker','o',...
    'LineStyle','none',...
    'Tag','strainlineInitialCondition')
if ~isfield(showPlot,'strainlineInitialCondition') || ...
        showPlot.strainlineInitialCondition == false
    set(findobj(axes,'tag','strainlineInitialCondition'),'visible',...
        'off')
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
    plot_strainline_segment(axes,strainline.position,...
        strainline.segmentIndex)
    if ~isfield(showPlot,'strainlineSegment') || ...
            showPlot.strainlineSegment == false
        set(findobj(axes,'tag','strainlineSegment'),'visible','off')
    end
end

if all(isfield(strainline,{'position','segmentIndex',...
        'filteredSegmentIndex'}))
    if showPlot.strainlineFiltered || showPlot.strainlineSegment
        
        switch strainline.filteringMethod
            case 'hausdorff'
                
            case 'superminimization'
                
                if showPlot.superminLine == true
                    showPlot.superminLine = axes;
                end
                
                if ~isfield(strainline,'filteredStrainlineIndex')
                    strainline.filteredStrainlineIndex = superminimize_grid(...
                        strainline.position,strainline.segmentIndex,...
                        strainline.relativeStretching,...
                        strainline.superminimizationDistance,...
                        flow.domain,strainline.resolution,showPlot.superminLine);
                end
                
            otherwise
                error('Filtering method not recognized')
        end
        
        if showPlot.strainlineFiltered
            plot_filtered_strainline(axes,strainline.position,...
                strainline.segmentIndex,strainline.filteredSegmentIndex)
        end
        
    end
end

function plot_geodesic_deviation_points(position,index,axes)
    
plot(axes,position(index,1),position(index,2),...
    'MarkerEdgeColor','r',...
    'Marker','o',...
    'LineStyle','none',...
    'Tag','geodesicDeviationPoint')

function plot_strainline(axes,strainline)

cellfun(@(position)plot(axes,position(:,1),position(:,2),'color','k',...
    'Tag','strainline'),strainline.position)

