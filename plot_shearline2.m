function hPlot = plot_shearline2(axes,position,domain,periodicBc)

if periodicBc(1)
    [positionPer,idx] = apply_periodic_bc(position,periodicBc,domain);
    if isempty(idx)
        hPlot = plot(axes,position(:,1),position(:,2));
    else
        idx = [0;idx;size(position,1)];
        nSegments = numel(idx)-1;
        hPlot = nan(nSegments,1);
        for m = 1:numel(idx)-1
            iIdx = idx(m)+1:idx(m+1);
            hPlot(m) = plot(axes,positionPer(iIdx,1),positionPer(iIdx,2));
        end
    end
else
    hPlot = plot(axes,position(:,1),position(:,2));
end
