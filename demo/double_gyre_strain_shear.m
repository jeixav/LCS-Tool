% Superimpose strain and shear structures
function hAxes = double_gyre_strain_shear

[doubleGyre,hStrainAxes] = double_gyre_strain;
delete(get(hStrainAxes,'parent'))

[~,hShearAxes,hDgAxes,averageGeodesicDeviationPos,closedOrbitPositionPos,averageGeodesicDeviationNeg,closedOrbitPositionNeg] = double_gyre_shear_closed;
delete(get(hDgAxes,'parent'))
delete(get(hShearAxes,'parent'))

hAxes = setup_figure(doubleGyre.flow.domain);
hPlot = plot_filtered_strainline(hAxes,doubleGyre.strainline.position,doubleGyre.strainline.segmentIndex,doubleGyre.strainline.filteredSegmentIndex);
set(hPlot,'color','k')
set(hPlot,'linewidth',1)

[~,idx] = min(averageGeodesicDeviationPos);
hPlot = plot(hAxes,closedOrbitPositionPos{idx}(:,1),closedOrbitPositionPos{idx}(:,2));
set(hPlot,'color','k')
set(hPlot,'linewidth',1)

[~,idx] = min(averageGeodesicDeviationNeg);
hPlot = plot(hAxes,closedOrbitPositionNeg{idx}(:,1),closedOrbitPositionNeg{idx}(:,2));
set(hPlot,'color','k')
set(hPlot,'linewidth',1)
