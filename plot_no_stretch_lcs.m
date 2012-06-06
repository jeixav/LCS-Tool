function plot_no_stretch_lcs(axes,flow,noStretchLine,showPlot)

if all(isfield(noStretchLine,{'chiPos','chiNeg'}))

    if ~isfield(showPlot,'etaQuiverDownsampling')
        showPlot.chiQuiverDownsampling = uint64(1);
    end

    plot_chi_quiver(axes,flow,noStretchLine,showPlot.chiQuiverDownsampling)
    
    if ~isfield(showPlot,'chiPosQuiver') || showPlot.chiPosQuiver == false
        set(findobj(axes,'tag','chiPosQuiver'),'visible','off')
    end
    
    if ~isfield(showPlot,'chiNegQuiver') || showPlot.chiNegQuiver == false
        set(findobj(axes,'tag','chiNegQuiver'),'visible','off')
    end

end

if all(isfield(noStretchLine,{'positionPos','positionNeg'}))
    
    plot_no_stretch_line(axes,noStretchLine)
    
    if ~isfield(showPlot,'noStretchLinePos') || ...
            showPlot.noStretchLinePos == false
        set(findobj(axes,'tag','noStretchLinePos'),'visible','off')
    end
    
    if ~isfield(showPlot,'noStretchLineNeg') || ...
            showPlot.noStretchLineNeg == false
        set(findobj(axes,'tag','noStretchLineNeg'),'visible','off')
    end
    
end

function plot_no_stretch_line(axes,noStretchLine)

cellfun(@(position)plot(axes,position(:,1),position(:,2),...
    'color','r','tag','noStretchLinePos'),...
    noStretchLine.positionPos);
cellfun(@(position)plot(axes,position(:,1),position(:,2),...
    'color','k','tag','noStretchLineNeg'),...
    noStretchLine.positionNeg);

function plot_chi_quiver(axes,flow,noStretchLine,downsampling)

dsResolution = flow.resolution/downsampling;
quiverScale = .5;

chiPosXGrid = reshape(noStretchLine.chiPos(:,1),fliplr(flow.resolution));
chiPosYGrid = reshape(noStretchLine.chiPos(:,2),fliplr(flow.resolution));
quiver(axes,...
    linspace(flow.domain(1,1),flow.domain(1,2),dsResolution(1)),...
    linspace(flow.domain(2,1),flow.domain(2,2),dsResolution(2)),...
    chiPosXGrid(1:downsampling:end,1:downsampling:end),...
    chiPosYGrid(1:downsampling:end,1:downsampling:end),...
    quiverScale,'color','r','tag','chiPosQuiver')

chiNegXGrid = reshape(noStretchLine.chiNeg(:,1),fliplr(flow.resolution));
chiNegYGrid = reshape(noStretchLine.chiNeg(:,2),fliplr(flow.resolution));
quiver(axes,...
    linspace(flow.domain(1,1),flow.domain(1,2),dsResolution(1)),...
    linspace(flow.domain(2,1),flow.domain(2,2),dsResolution(2)),...
    chiNegXGrid(1:downsampling:end,1:downsampling:end),...
    chiNegYGrid(1:downsampling:end,1:downsampling:end),...
    quiverScale,'color','k','tag','chiNegQuiver')
