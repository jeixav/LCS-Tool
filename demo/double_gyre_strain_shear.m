% Superimpose strain and shear structures
function hAxes = double_gyre_strain_shear

[doubleGyre,hStrainAxes] = double_gyre_strain;
delete(get(hStrainAxes,'parent'))

[~,hShearAxes,closedOrbitPositionPos,closedOrbitPositionNeg] = double_gyre_shear_closed;
delete(get(hShearAxes,'parent'))

hAxes = setup_figure(doubleGyre.flow.domain);
hPlot = plot_filtered_strainline(hAxes,doubleGyre.strainline.position,doubleGyre.strainline.segmentIndex,doubleGyre.strainline.filteredSegmentIndex);
set(hPlot,'color','r')
set(hPlot,'linewidth',1)

hPlot = plot(hAxes,closedOrbitPositionPos{end}(:,1),closedOrbitPositionPos{end}(:,2));
set(hPlot,'color','g')
set(hPlot,'linewidth',1)

hPlot = plot(hAxes,closedOrbitPositionNeg{1}(:,1),closedOrbitPositionNeg{1}(:,2));
set(hPlot,'color','g')
set(hPlot,'linewidth',1)
