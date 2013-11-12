% plot_ftle Plot finite-time Lyapunov exponent
%
% SYNTAX
% [hFtle,hColorbar] = plot_ftle(hAxes,flow,ftle)
%
% DESCRIPTION
% To adjust FTLE range: set(hAxes,'clim',[0,.5]);

function [hFtle,hColorbar] = plot_ftle(hAxes,domain,resolution,ftle)

ftle = reshape(ftle,fliplr(resolution));

hFtle = imagesc(domain(1,:),domain(2,:),ftle);
set(hFtle,'Parent',hAxes);
set(hFtle,'tag','ftle')

set(hAxes,'layer','top')

hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','FTLE')
