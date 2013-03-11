function plot_ftle(hAxes,flow,ftle)

ftle = reshape(ftle,fliplr(flow.resolution));

% nContour = 40;
% contourMin = min(ftle(:));
% contourMax = max(ftle(:));

h = imagesc(flow.domain(1,:),flow.domain(2,:),ftle);
set(h,'Parent',hAxes);
set(h,'tag','ftle')

colorbar
