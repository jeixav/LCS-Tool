%% Demo script of closed orbit detection
% Demo of poincare_closed_orbit_multi()
% Detection of the outermost closed shearline for multiple poincare
% sections
%
% fhuhn - 2013/07/11
% last changed: 2013/07/11

lcs_toolbox_path = fullfile('..','..');

% add toolbox folder
addpath(lcs_toolbox_path);
addpath(fullfile(lcs_toolbox_path,'ParforProgress2'));

% open matlab pool for parallel computing
if ~matlabpool('size')
    matlabpool('open');
end

% FIXME Provide commands to compute dataset
load run02_ocean_150x150

%% plot FTLE for orientation
lambda2 = reshape(oceanAgulhas.flow.cgEigenvalue(:,2), fliplr(oceanAgulhas.flow.resolution));

hFigure = figure;
hAxes = axes('parent',hFigure);
set(hAxes,'nextplot','add')
hImagesc = imagesc(lon_axis, lat_axis, log10(sqrt(lambda2)));
set(hImagesc,'parent',hAxes)
axis(hAxes,'equal')
axis(hAxes,'tight')
xlabel(hAxes,'Lon')
ylabel(hAxes,'Lat')

%% define poincare sections
% Poincare section should be placed with 1st point in center of elliptic
% region and with second point outside the elliptic region
%
% poincareSection{i}.endPosition = [x1 y1; x2 y2]
% poincareSection{i}.numPoints = nPoints
% poincareSection{i}.integrationLength = integrationLength

% Example with 2 poincare sections
% poincare section 1
poincareSection{1}.endPosition = [-15.4661 -29.3477; -15.6653 -29.9426];
% poincare section 2
poincareSection{2}.endPosition = [-17 -32; -16.5 -31.5];

nPoincareSection = numel(poincareSection);
rOrbits = nan(1,nPoincareSection);
integrationLength = nan(1,nPoincareSection);
for i=1:nPoincareSection
    poincareSection{i}.numPoints = 400; %#ok<SAGROW>
    % radius = length of poincare section
    rOrbits(i) = sqrt((poincareSection{i}.endPosition(2,1)-poincareSection{i}.endPosition(1,1)).^2 + ...
        (poincareSection{i}.endPosition(2,2)-poincareSection{i}.endPosition(1,2)).^2 );
    % set integration length conservatively = twice the expected circumference
    poincareSection{i}.integrationLength = [0 2*(2*pi*rOrbits(i))]; %#ok<SAGROW>
end

arrayfun(@(idx)plot(hAxes,poincareSection{idx}.endPosition(:,1), poincareSection{idx}.endPosition(:,2), 'ko-', 'linewidth',2), 1:size(poincareSection,2) );
drawnow

%% closed orbit detection
odeSolverOptions = odeset('relTol',1e-6);
showPlot = 1;
[closedOrbits, orbits] = poincare_closed_orbit_multi(oceanAgulhas.flow,oceanAgulhas.shearline,poincareSection,odeSolverOptions,showPlot);

%% plot orbits
hFigure = figure;
hAxes = axes;
set(hAxes,'nextplot','add')
hImagesc = imagesc(lon_axis, lat_axis, log10(sqrt(lambda2)));
set(hImagesc,'parent',hAxes)
% poincare sections
arrayfun(@(i)plot(hAxes,poincareSection{i}.endPosition(:,1),poincareSection{i}.endPosition(:,2),'ro--','linewidth',2),1:nPoincareSection);
% orbits
arrayfun(@(i)plot(hAxes,orbits{1}{1}{i}(:,1),orbits{1}{1}{i}(:,2),'b-'),1:poincareSection{1}.numPoints);
arrayfun(@(i)plot(hAxes,orbits{1}{2}{i}(:,1),orbits{1}{2}{i}(:,2),'b-'),1:poincareSection{1}.numPoints);
arrayfun(@(i)plot(hAxes,orbits{2}{1}{i}(:,1),orbits{2}{1}{i}(:,2),'b-'),1:poincareSection{2}.numPoints);
arrayfun(@(i)plot(hAxes,orbits{2}{2}{i}(:,1),orbits{2}{2}{i}(:,2),'b-'),1:poincareSection{2}.numPoints);
% etaPos closed orbits
arrayfun(@(i)plot(hAxes,closedOrbits{i}{1}(:,1),closedOrbits{i}{1}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
% etaNeg closed orbits
arrayfun(@(i)plot(hAxes,closedOrbits{i}{2}(:,1),closedOrbits{i}{2}(:,2),'color',[0 0.6 0],'linewidth',2),1:size(closedOrbits,2));
colormap(gray)
colorbar('peer',hAxes)
axis(hAxes,'equal')
axis(hAxes,'tight')
