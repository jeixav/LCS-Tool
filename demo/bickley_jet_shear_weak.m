function bickleyJet = bickley_jet_shear_weak

bickleyJet = bickley_jet(3);

bickleyJet.flow = set_flow_resolution(1000,bickleyJet.flow);
bickleyJet.flow = set_flow_ode_solver_options(odeset('RelTol',1e-4),bickleyJet.flow);
bickleyJet.flow.cgStrainMethod.name = 'equationOfVariation';
bickleyJet.flow.cgStrainMethod.eigenvalueFromMainGrid = false;
bickleyJet.flow.cgStrainMethod.auxiliaryGridRelativeDelta = 1e-2;
bickleyJet.flow.cgStrainCustomEigMethod = false;
bickleyJet.flow.imposeIncompressibility = true;

load('datasets/bickley_jet/bickleyJet3_t20_2.mat')

bickleyJet.shearline = set_shearline_resolution([20,10],bickleyJet.shearline);
bickleyJet.shearline = set_shearline_max_length(2e7,bickleyJet.shearline);
bickleyJet.shearline = set_shearline_ode_solver_options(odeset('RelTol',1e-4),bickleyJet.shearline);
bickleyJet.shearline = set_shearline_average_geodesic_deviation_tol([1,1]*5e-5,bickleyJet.shearline);

showPlot.shearlinePos = true;
showPlot.shearlinePosFiltered = true;
showPlot.shearlineNegFiltered = false;

verbose.progress = false;

[bickleyJet,hAxes] = shear_lcs_script(bickleyJet,showPlot,verbose);
set(findobj(hAxes,'tag','shearlinePosFiltered'),'visible','off')

%% Poincare analysis
poincareSection.endPosition = [4.5e6,-1.5e6;6e6,-1.5e6];
plot_poincare_section(hAxes,poincareSection.endPosition);
poincareSection.numPoints = 400;
poincareSection.tmax = 2e7;
[~,closedOrbitPositionPos] = poincare_closed_orbit(bickleyJet.flow,bickleyJet.shearline.etaPos,poincareSection,bickleyJet.shearline.odeSolverOptions,hAxes);
set(findobj(hAxes,'tag','closedOrbit'),'visible','off')

%% Highlight outermost closed orbit
hClosedOutermost = plot(hAxes,closedOrbitPositionPos{1}(:,1),closedOrbitPositionPos{1}(:,2));
set(hClosedOutermost,'color','g')
set(hClosedOutermost,'linewidth',1)
set(hClosedOutermost,'tag','closedOrbitOutermost')

% plot_negative_shearlines(bickleyJet)

function plot_negative_shearlines(bickleyJet)

showPlot.shearlineNeg = true;
showPlot.shearlinePosFiltered = false;
showPlot.shearlineNegFiltered = true;

shear_lcs_script(bickleyJet,showPlot)

function hPoincareSection = plot_poincare_section(hAxes,poincareSectionPosition)

hPoincareSection = plot(hAxes,poincareSectionPosition(:,1),poincareSectionPosition(:,2));

set(hPoincareSection,'tag','poincareSection')
set(hPoincareSection,'color','k')
set(hPoincareSection,'linestyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerSize',2)
set(hPoincareSection,'MarkerFaceColor','k')
set(hPoincareSection,'MarkerEdgeColor','k')
