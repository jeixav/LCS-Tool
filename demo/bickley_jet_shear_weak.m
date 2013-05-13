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
showPlot.shearlinePosFiltered = false;
showPlot.shearlineNegFiltered = false;

verbose.progress = false;

load('datasets/bickley_jet/bickleyJet3_t20_3.mat')
[bickleyJet,hAxes] = shear_lcs_script(bickleyJet,showPlot,verbose);

%% Manually highlight closed orbits
h = plot(hAxes,bickleyJet.shearline.positionPos{52}(:,1),bickleyJet.shearline.positionPos{52}(:,2));
set(h,'color','g')
set(h,'linewidth',1)

h = plot(hAxes,bickleyJet.shearline.positionPos{114}(:,1),bickleyJet.shearline.positionPos{114}(:,2));
set(h,'color','g')
set(h,'linewidth',1)

h = plot(hAxes,bickleyJet.shearline.positionPos{20}(:,1),bickleyJet.shearline.positionPos{20}(:,2));
set(h,'color','g')
set(h,'linewidth',1)
set(h,'linestyle','--')

h = plot(hAxes,bickleyJet.shearline.positionPos{67}(:,1),bickleyJet.shearline.positionPos{67}(:,2));
set(h,'color','g')
set(h,'linewidth',1)
set(h,'linestyle','--')

h = plot(hAxes,bickleyJet.shearline.positionPos{2}(:,1),bickleyJet.shearline.positionPos{2}(:,2));
set(h,'color','g')
set(h,'linewidth',1)

h = plot(hAxes,bickleyJet.shearline.positionPos{192}(:,1),bickleyJet.shearline.positionPos{192}(:,2));
set(h,'color','g')
set(h,'linewidth',1)

% %% Poincare analysis
% poincareSection.endPosition = [4.5e6,-1.5e6;6e6,-1.5e6];
% plot_poincare_section(hAxes,poincareSection.endPosition);
% poincareSection.numPoints = 400;
% poincareSection.tmax = 2e7;
% [~,closedOrbitPositionPos] = poincare_closed_orbit(bickleyJet.flow,bickleyJet.shearline.etaPos,poincareSection,bickleyJet.shearline.odeSolverOptions,hAxes);
% set(findobj(hAxes,'tag','closedOrbit'),'visible','off')
% 
% %% Highlight outermost closed orbit
% hClosedOutermost = plot(hAxes,closedOrbitPositionPos{1}(:,1),closedOrbitPositionPos{1}(:,2));
% set(hClosedOutermost,'color','g')
% set(hClosedOutermost,'linewidth',1)
% set(hClosedOutermost,'tag','closedOrbitOutermost')

[bickleyJet,hAxes] = plot_negative_shearlines(bickleyJet);

jetLowerIdx = 11; % 32 51 81 85 92 111 141 162 181
hJetLower = plot(hAxes,bickleyJet.shearline.positionNeg{jetLowerIdx}(:,1),bickleyJet.shearline.positionNeg{jetLowerIdx}(:,2));
set(hJetLower,'color','g')
set(hJetLower,'linewidth',1)
set(hJetLower,'linestyle','--')

closed1Idx = 17; % 49 46 29
hClosed1 = plot(hAxes,bickleyJet.shearline.positionNeg{closed1Idx}(:,1),bickleyJet.shearline.positionNeg{closed1Idx}(:,2));
set(hClosed1,'color','g')
set(hClosed1,'linewidth',1)

jetCentreIdx = 34; % 56 94 104 115 146 164
hJetCentre = plot(hAxes,bickleyJet.shearline.positionNeg{jetCentreIdx}(:,1),bickleyJet.shearline.positionNeg{jetCentreIdx}(:,2));
set(hJetCentre,'color','g')
set(hJetCentre,'linewidth',1)
set(hJetCentre,'linestyle','--')

closed2Idx = 119; % 89 117 118 99 88 96 97 98 106 107 108 87 109
hClosed2 = plot(hAxes,bickleyJet.shearline.positionNeg{closed2Idx}(:,1),bickleyJet.shearline.positionNeg{closed2Idx}(:,2));
set(hClosed2,'color','g')
set(hClosed2,'linewidth',1)

closed3Idx = 156; % 187 188 159 169 179
hClosed3 = plot(hAxes,bickleyJet.shearline.positionNeg{closed3Idx}(:,1),bickleyJet.shearline.positionNeg{closed3Idx}(:,2));
set(hClosed3,'color','g')
set(hClosed3,'linewidth',1)

function [bickleyJet,hAxes] = plot_negative_shearlines(bickleyJet)

showPlot.shearlineNeg = true;
showPlot.shearlinePosFiltered = false;
showPlot.shearlineNegFiltered = false;

[bickleyJet,hAxes] = shear_lcs_script(bickleyJet,showPlot);

function hPoincareSection = plot_poincare_section(hAxes,poincareSectionPosition)

hPoincareSection = plot(hAxes,poincareSectionPosition(:,1),poincareSectionPosition(:,2));

set(hPoincareSection,'tag','poincareSection')
set(hPoincareSection,'color','k')
set(hPoincareSection,'linestyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerSize',2)
set(hPoincareSection,'MarkerFaceColor','k')
set(hPoincareSection,'MarkerEdgeColor','k')
