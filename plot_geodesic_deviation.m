function hAxes = plot_geodesic_deviation(geodesicDeviation,domain)

hFigure = figure;
hAxes = axes;
set(hAxes,'parent',hFigure)
set(hAxes,'nextplot','add')
set(hAxes,'DataAspectRatioMode','manual')
set(hAxes,'DataAspectRatio',[1,1,1])
set(hAxes,'xlim',domain(1,:))
set(hAxes,'ylim',domain(2,:))
hXlabel = xlabel('x');
set(hXlabel,'parent',hAxes)
hYlabel = ylabel('y');
set(hYlabel,'parent',hAxes)
hImagesc = imagesc(domain(1,:),domain(2,:),log10(geodesicDeviation));
set(hImagesc,'parent',hAxes);
uistack(hImagesc,'bottom')
cbar_axes = colorbar;
set(cbar_axes,'parent',get(hAxes,'parent'))
set(get(cbar_axes,'xlabel'),'string','log(d_g)')
