% plot_shear_lcs Plot shearline LCSs: elliptic and parabolic barriers.
%
% DESCRIPTION
% Every object in the plot has a tag. This can be used to control
% visibility and to set properties.
% unique(get(get(gca,'children'),'tag'))
% lists all object tags.
%
% set(findobj(gca,'tag','shearlinePosFiltered'),'visible','off')
% Hides filtered eta+ shearlines.

function plot_shear_lcs(axes,flow,shearline,showPlot)

if all(isfield(shearline,{'etaPos','etaNeg'}))

    if ~isfield(showPlot,'etaQuiverDownsampling')
        showPlot.etaQuiverDownsampling = uint64(1);
    end
    
    plot_eta_quiver(axes,flow,shearline,showPlot.etaQuiverDownsampling)
    
    if ~isfield(showPlot,'etaPosQuiver') || showPlot.etaPosQuiver == false
        set(findobj(axes,'tag','etaPosQuiver'),'visible','off')
    end
    
    if ~isfield(showPlot,'etaNegQuiver') || showPlot.etaNegQuiver == false
        set(findobj(axes,'tag','etaNegQuiver'),'visible','off')
    end
    
end

if all(isfield(shearline,{'positionPos','positionNeg'}))
    
    [hShearlinePos,hShearlineNeg] = plot_shearline(axes,flow,shearline);
    set([hShearlinePos;hShearlineNeg],'color',.8*[1 1 1])

    if ~isfield(showPlot,'shearlinePos') ...
            || showPlot.shearlinePos == false
        set(hShearlinePos,'visible','off')
    end
    
    if ~isfield(showPlot,'shearlineNeg') ...
            || showPlot.shearlineNeg == false
        set(hShearlineNeg,'visible','off')
    end
    
end

if all(isfield(shearline,{'positionClosedPos','positionClosedNeg'}))
    
    plot_closed_shearline(axes,shearline)
    
    if ~isfield(showPlot,'shearlinePosClosed') ...
            || showPlot.shearlinePosClosed == false
        set(findobj(axes,'tag','shearlinePosClosed'),'visible','off')
    end
    
    if ~isfield(showPlot,'shearlineNegClosed') ...
            || showPlot.shearlineNegClosed == false
        set(findobj(axes,'tag','shearlineNegClosed'),'visible','off')
    end

end

if all(isfield(shearline,{'filteredIndexPos','filteredIndexPos'}))

    plot_filtered_shearline(axes,shearline)

    if ~isfield(showPlot,'shearlinePosFiltered') ...
            || showPlot.shearlinePosFiltered == false
        set(findobj(axes,'tag','shearlinePosFiltered'),'visible','off')
    end
    
    if ~isfield(showPlot,'shearlineNegFiltered') ...
            || showPlot.shearlineNegFiltered == false
        set(findobj(axes,'tag','shearlineNegFiltered'),'visible','off')
    end
    
end

if isfield(shearline,'initialPosition')
    
    plot_shearline_initial_position(axes,shearline)
    
    if ~isfield(showPlot,'shearlineInitialPosition') || ...
            showPlot.shearlineInitialPosition == false
        set(findobj(axes,'tag','shearlineInitialPosition'),'visible','off')
    end
    
end

function plot_shearline_initial_position(axes,shearline)

hPlot = plot(axes,shearline.initialPosition(:,1),shearline.initialPosition(:,2));

set(hPlot,'tag','shearlineInitialPosition')
set(hPlot,'MarkerFaceColor','k')
set(hPlot,'MarkerEdgeColor','k')
set(hPlot,'Marker','o')
set(hPlot,'LineStyle','none')

function plot_filtered_shearline(axes,shearline)

cellfun(@(position)plot(axes,position(:,1),position(:,2),...
    'color','r','tag','shearlinePosFiltered'),...
    shearline.positionPos(shearline.filteredIndexPos));
cellfun(@(position)plot(axes,position(:,1),position(:,2),...
    'color','k','tag','shearlineNegFiltered'),...
    shearline.positionNeg(shearline.filteredIndexNeg));

function [hShearlinePos,hShearlineNeg] = plot_shearline(axes,flow,shearline)

hShearlinePos = cellfun(@(position)plot_shearline2(axes,position,flow.domain,flow.periodicBc),shearline.positionPos,'UniformOutput',false);
hShearlinePos = vertcat(hShearlinePos{:});
set(hShearlinePos,'tag','shearlinePos')
set(hShearlinePos,'color','r')

hShearlineNeg = cellfun(@(position)plot_shearline2(axes,position,flow.domain,flow.periodicBc),shearline.positionNeg,'UniformOutput',false);
hShearlineNeg = vertcat(hShearlineNeg{:});
set(hShearlineNeg,'tag','shearlineNeg')
set(hShearlineNeg,'color','k')

function plot_closed_shearline(axes,shearline)

cellfun(@(position)plot(axes,position(:,1),position(:,2),...
    'color','r','linestyle','--','tag','shearlinePosClosed'),...
    shearline.positionClosedPos);
cellfun(@(position)plot(axes,position(:,1),position(:,2),...
    'color','k','linestyle','--','tag','shearlineNegClosed'),...
    shearline.positionClosedNeg);

% FIXME plot_strain_lcs does not use a separate function to produce similar
% plot.
function plot_eta_quiver(hAxes,flow,shearline,downsampling)

dsResolution = flow.resolution/downsampling;
quiverScale = .5;

etaPosXGrid = reshape(shearline.etaPos(:,1),fliplr(flow.resolution));
etaPosYGrid = reshape(shearline.etaPos(:,2),fliplr(flow.resolution));
quiver(hAxes,...
    linspace(flow.domain(1,1),flow.domain(1,2),dsResolution(1)),...
    linspace(flow.domain(2,1),flow.domain(2,2),dsResolution(2)),...
    etaPosXGrid(1:downsampling:end,1:downsampling:end),...
    etaPosYGrid(1:downsampling:end,1:downsampling:end),...
    quiverScale,'color','r','tag','etaPosQuiver')

etaNegXGrid = reshape(shearline.etaNeg(:,1),fliplr(flow.resolution));
etaNegYGrid = reshape(shearline.etaNeg(:,2),fliplr(flow.resolution));
quiver(hAxes,...
    linspace(flow.domain(1,1),flow.domain(1,2),dsResolution(1)),...
    linspace(flow.domain(2,1),flow.domain(2,2),dsResolution(2)),...
    etaNegXGrid(1:downsampling:end,1:downsampling:end),...
    etaNegYGrid(1:downsampling:end,1:downsampling:end),...
    quiverScale,'color','k','tag','etaNegQuiver')
