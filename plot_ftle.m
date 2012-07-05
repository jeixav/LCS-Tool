function plot_ftle(axes,flow,ftle)

ftle = reshape(ftle,fliplr(flow.resolution));

nContour = 40;
contourMin = min(ftle(:));
contourMax = max(ftle(:));

[~,h] = contourf(axes,linspace(flow.domain(1,1),flow.domain(1,2),...
    flow.resolution(1)),linspace(flow.domain(2,1),flow.domain(2,2),...
    flow.resolution(2)),ftle,linspace(contourMin,contourMax,nContour),...
    'linestyle','none');
set(h,'tag','ftle')

colorbar
