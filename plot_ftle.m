function [hFtle,hColorbar] = plot_ftle(hAxes,flow,ftle)

ftle = reshape(ftle,fliplr(flow.resolution));

hFtle = imagesc(flow.domain(1,:),flow.domain(2,:),ftle);
set(hFtle,'Parent',hAxes);
set(hFtle,'tag','ftle')

hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
