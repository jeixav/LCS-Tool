% poincare_closed_orbit_mod(flow,vectorField,poincareSection,...
% odeSolverOptions,timespan,showGraph)
%
% Find closed orbit using Poincare section
%
% EXAMPLE
% addpath('flow_templates')
% matlabpool('open')
% pctRunOnAll javaaddpath('ParforProgress2')
%
% bickleyJet = bickley_jet(1);
%
% bickleyJet.flow = set_flow_resolution(40,bickleyJet.flow);
% bickleyJet.flow = set_flow_ode_solver_options(odeset('relTol',1e-4,...
%     'absTol',1e-6),bickleyJet.flow);
% bickleyJet.flow = set_flow_domain([3.5 4.7; -1.3 -.2],bickleyJet.flow);
% bickleyJet.flow.imposeIncompressibility = false;
%
% bickleyJet.shearline = rmfield(bickleyJet.shearline,'odeSolverOptions');
% bickleyJet.shearline = set_shearline_resolution([4 4],...
%     bickleyJet.shearline);
%
% showPlot.shearlinePosFiltered = false;
% showPlot.shearlineNegFiltered = false;
% showPlot.etaPosQuiver = true;
%
% bickleyJet = shear_lcs_script(bickleyJet,showPlot);
%
% poincareSection.endPosition = [4.15 -.6; 4.15 -.3];
% poincareSection.numPoints = 50;
% odeSolverOptions = odeset('relTol',1e-6);
%
% closedOrbitInitialPosition = poincare_closed_orbit(bickleyJet.flow,...
%     bickleyJet.shearline.etaPos,poincareSection,odeSolverOptions);
%
% disp('Closed orbit positions:')
% disp(num2str(transpose(closedOrbitInitialPosition)))

function [closedOrbitPosition, orbitPosition] = poincare_closed_orbit_mod(flow,...
    vectorField,poincareSection,odeSolverOptions,timespan,showGraph)

narginchk(5,6)

if nargin == 5
    showGraph = true;
end

% Initial positions for Poincare orbits
orbitInitialPositionX = linspace(poincareSection.endPosition(1,1),...
    poincareSection.endPosition(2,1),poincareSection.numPoints);
orbitInitialPositionY = linspace(poincareSection.endPosition(1,2),...
    poincareSection.endPosition(2,2),poincareSection.numPoints);
orbitInitialPosition = transpose([orbitInitialPositionX; ...
    orbitInitialPositionY]);

% plot all initial positions
% if showGraph
%     hPoincareSection = plot(orbitInitialPosition(:,1),...
%         orbitInitialPosition(:,2),'-x', 'tag', 'orbitInitialPosition');
%     hParent = get(hPoincareSection,'parent');
% end

flowDomain = flow.domain;
flowResolution = flow.resolution;
orbitPosition = cell(poincareSection.numPoints,1);

% integrate orbits
parfor idx = 1:poincareSection.numPoints
    orbitPosition{idx} = integrate_line_closed_mod(timespan,...
        orbitInitialPosition(idx,:),flowDomain,flowResolution,...
        vectorField,poincareSection.endPosition,odeSolverOptions);
end

% plot all orbits
% if showGraph
%     arrayfun(@(idx)plot(hParent,orbitPosition{idx}(:,1),...
%         orbitPosition{idx}(:,2), 'tag', 'orbitPosition'),1:poincareSection.numPoints);
% end

% final position of orbits
orbitFinalPosition = cellfun(@(position)position(end,:),orbitPosition,...
    'UniformOutput',false);
orbitFinalPosition = cell2mat(orbitFinalPosition);

if showGraph
    hfigure = figure;
    hAxes = axes;
    set(hAxes,'parent',hfigure);
    set(hAxes,'nextplot','add')
    set(hAxes,'box','on');
    set(hAxes,'xgrid','on');
    set(hAxes,'ygrid','on');    
    xLength = sqrt(diff(poincareSection.endPosition(:,1))^2 ...
        + diff(poincareSection.endPosition(:,2))^2);
end

% angle of poincare section with x-axis
theta = -atan((poincareSection.endPosition(1,2) - poincareSection.endPosition(2,2))...
    /(poincareSection.endPosition(1,1) - poincareSection.endPosition(2,1)));
rotationMatrix = [cos(theta),-sin(theta);sin(theta),cos(theta)];

% Translate to origin
s(:,1) = orbitInitialPosition(:,1) - poincareSection.endPosition(1,1);
s(:,2) = orbitInitialPosition(:,2) - poincareSection.endPosition(1,2);

t(:,1) = orbitFinalPosition(:,1) - poincareSection.endPosition(1,1);
t(:,2) = orbitFinalPosition(:,2) - poincareSection.endPosition(1,2);

% Rotate to 0
s = rotationMatrix*transpose(s);
t = rotationMatrix*transpose(t);
s = transpose(s);
t = transpose(t);

if showGraph
    set(hAxes,'xlim',[0 xLength]);
    plot(hAxes,abs(s(:,1)),t(:,1)-s(:,1),'-x');
    title('Closed orbit detection - Poincare section - return distance')
    xlabel(hAxes,'s');
    ylabel(hAxes,'p(s) - s');
end


% find zero crossings of poincare return map (linear)
[~,closedOrbitInitialPosition] = crossing(t(:,1) - s(:,1),s(:,1));


if isempty(closedOrbitInitialPosition)
    
    closedOrbitInitialPosition(1,1:2) = [NaN NaN];
    closedOrbitPosition = [NaN NaN];
    
else
    
    if showGraph
        [nClosedOrbit ~] = size(closedOrbitInitialPosition);
        h1 = plot(hAxes,abs(closedOrbitInitialPosition),zeros(1,nClosedOrbit),'ro');
        legend([h1(1) 0], 'Closed orbits', 'No valid closed orbits');
    end
    
    % Rotate to theta
    xx = [transpose(closedOrbitInitialPosition) ...
        zeros(numel(closedOrbitInitialPosition),1)];
    xx = rotationMatrix\transpose(xx);
    xx = transpose(xx);
    
    % Translate from origin
    closedOrbitInitialPositionX = xx(:,1) + poincareSection.endPosition(1,1);
    closedOrbitInitialPositionY = xx(:,2) + poincareSection.endPosition(1,2);
    
    closedOrbitInitialPosition = [closedOrbitInitialPositionX,...
        closedOrbitInitialPositionY];
        
    % FILTER: Discard closed orbits if images of neighbor points do not fall on poincare section
    % i.e. discard zero crossings due to outlyers
    %***********************
    alphaThresh = 1e-3;
    %***********************
    [nClosedOrbit ~] = size(closedOrbitInitialPosition);
    for i=1:nClosedOrbit
        % find 2 neighbor points
        [orbitInitialPositionSorted, ix] = sort(orbitInitialPosition(:,1));
        indx10 = max(find( closedOrbitInitialPosition(i,1) > orbitInitialPositionSorted ));
        indx20 = min(find( closedOrbitInitialPosition(i,1) < orbitInitialPositionSorted ));
        indx1 = min( ix(indx10), ix(indx20));
        indx2 = max( ix(indx10), ix(indx20));
        if indx2 <= indx1 || abs(indx1-indx2) ~=1
            error('Selection of neighbor orbits failed.')
        end
        % check if image of neighbor points falls on poincare section
        % i.e. angle of vectors to endpoints of poincare section = pi
        % lower neighbor
        v1 = poincareSection.endPosition(1,:) - orbitFinalPosition(indx1,:);
        v2 = poincareSection.endPosition(2,:) - orbitFinalPosition(indx1,:);
        alpha1 = angle2vectors(v1, v2);
        % upper neighbor
        v1 = poincareSection.endPosition(1,:) - orbitFinalPosition(indx2,:);
        v2 = poincareSection.endPosition(2,:) - orbitFinalPosition(indx2,:);
        alpha2 = angle2vectors(v1, v2);
        % both angles have to be very close to pi
        if abs(alpha1-pi) > alphaThresh || abs(alpha2-pi) > alphaThresh
            closedOrbitInitialPosition(i,:) = NaN;
        end
    end
    % erase erroneous closed orbits
    [iy ~]= find(isnan(closedOrbitInitialPosition));
    closedOrbitInitialPosition(unique(iy),:) = [];
    
    
    if ~isempty(closedOrbitInitialPosition)
        [nClosedOrbit ~] = size(closedOrbitInitialPosition);
        % integrate closed orbits
        parfor idx = 1:nClosedOrbit
            closedOrbitPosition{idx} = integrate_line_closed_mod(timespan,...
                closedOrbitInitialPosition(idx,:),flowDomain,flowResolution,...
                vectorField,poincareSection.endPosition,odeSolverOptions);
        end
        
        %         if showGraph
        %             if ~isempty(closedOrbitInitialPositionY)
        %                 hClosedOrbit = arrayfun(@(idx)plot(hParent,...
        %                     closedOrbitPosition{idx}(:,1),...
        %                     closedOrbitPosition{idx}(:,2), 'tag','closedOrbitPosition'),1:nClosedOrbit);
        %                 set(hClosedOrbit,'color','r')
        %                 set(hClosedOrbit,'linewidth',2)
        %             end
        %         end
        
        % FILTER: select outermost closed orbit
        s1(:,1) = closedOrbitInitialPosition(:,1) - poincareSection.endPosition(1,1);
        s1(:,2) = closedOrbitInitialPosition(:,2) - poincareSection.endPosition(1,2);
        distR = sqrt(s1(:,1).^2 + s1(:,2).^2);
        % outermost = largest distance from 1st point of poincare section
        indR = find( distR == max(distR) );
        closedOrbitInitialPosition = closedOrbitInitialPosition(indR,:);
        closedOrbitPosition = closedOrbitPosition{indR}(:,:);
        
        if showGraph
            [nClosedOrbit ~] = size(closedOrbitInitialPosition);
            h2 = plot(hAxes,distR(indR),0,'go', 'markersize', 10);
            legend([h1(1) h2(1)], 'Closed orbits', 'Outermost valid closed orbit');
        end
        
    else
        closedOrbitInitialPosition(1,1:2) = [NaN NaN];
        closedOrbitPosition = [NaN NaN];
    end
    
end


function alpha = angle2vectors(v1, v2)
v1_norm = sqrt(v1(1)^2 + v1(2)^2);
v2_norm = sqrt(v2(1)^2 + v2(2)^2);
alpha = acos( (v1(1)*v2(1)+v1(2)*v2(2))/(v1_norm*v2_norm) );

