function [doubleGyre,hShearAxes,closedOrbitPositionPos,closedOrbitPositionNeg] = double_gyre_shear_closed

doubleGyre = double_gyre;

doubleGyre.flow = set_flow_resolution(500,doubleGyre.flow);
doubleGyre.flow = set_flow_timespan([0,20],doubleGyre.flow);
doubleGyre.flow = set_flow_ode_solver_options(odeset('RelTol',1e-8,'AbsTol',1e-8),doubleGyre.flow);
doubleGyre.flow.cgStrainMethod.name = 'finiteDifference';
doubleGyre.flow.coupledIntegration = 1e4;

doubleGyre.shearline = set_shearline_resolution([10,5],doubleGyre.shearline);
doubleGyre.shearline = set_shearline_max_length(4,doubleGyre.shearline);
doubleGyre.shearline = set_shearline_average_geodesic_deviation_tol([.32,.34],doubleGyre.shearline);

%% Positive shearlines
showPlot.shearlinePos = true;
showPlot.shearlinePosFiltered = false;
showPlot.shearlineNegFiltered = false;

verbose.progress = false;

[doubleGyre,hShearAxes] = shear_lcs_script(doubleGyre,showPlot,verbose);

poincareSection.endPosition = [.6,.55;.75,.55];
hPoincareSection = plot(hShearAxes,poincareSection.endPosition(:,1),poincareSection.endPosition(:,2));
set(hPoincareSection,'tag','poincareSection')
set(hPoincareSection,'color','k')
set(hPoincareSection,'linestyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerSize',2)
set(hPoincareSection,'MarkerFaceColor','k')
set(hPoincareSection,'MarkerEdgeColor','k')

poincareSection.numPoints = 200;
odeSolverOptions = odeset(doubleGyre.shearline.odeSolverOptions);

[~,closedOrbitPositionPos] = poincare_closed_orbit(doubleGyre.flow,doubleGyre.shearline.etaPos,poincareSection,odeSolverOptions,hShearAxes);
set(findobj(hShearAxes,'tag','closedOrbit'),'linewidth',.5)

%% Highlight longest closed orbit
hPlot = plot(hShearAxes,closedOrbitPositionPos{end}(:,1),closedOrbitPositionPos{end}(:,2));
set(hPlot,'tag','closedOrbitMinDg')
set(hPlot,'color','g')
set(hPlot,'linewidth',1)

%% Negative shearlines
clear('showPlot')
showPlot.shearlineNeg = true;
showPlot.shearlinePosFiltered = false;
showPlot.shearlineNegFiltered = false;

[doubleGyre,hShearAxes] = shear_lcs_script(doubleGyre,showPlot,verbose);

poincareSection.endPosition = [1.25,.45;1.45,.45];
hPoincareSection = plot(hShearAxes,poincareSection.endPosition(:,1),poincareSection.endPosition(:,2));
set(hPoincareSection,'tag','poincareSection')
set(hPoincareSection,'color','k')
set(hPoincareSection,'linestyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'markersize',2)
set(hPoincareSection,'MarkerFaceColor','k')
set(hPoincareSection,'MarkerEdgeColor','k')

poincareSection.numPoints = 200;
odeSolverOptions = odeset(doubleGyre.shearline.odeSolverOptions);

[~,closedOrbitPositionNeg,hPsAxes] = poincare_closed_orbit(doubleGyre.flow,doubleGyre.shearline.etaNeg,poincareSection,odeSolverOptions,hShearAxes);

set(findobj(hShearAxes,'tag','closedOrbit'),'linewidth',.5)

% Manually restrict Poincare return map axis limits
set(hPsAxes,'ylim',[-.0005,.0005])

%% Highlight longest closed orbit
hPlot = plot(hShearAxes,closedOrbitPositionNeg{1}(:,1),closedOrbitPositionNeg{1}(:,2));
set(hPlot,'tag','closedOrbitMinDg')
set(hPlot,'color','g')
set(hPlot,'linewidth',1)
