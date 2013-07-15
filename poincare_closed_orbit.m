% poincare_closed_orbit(flow,vectorField,poincareSection,...
% odeSolverOptions,timespan,showGraph)
%
% Find closed orbit using Poincare section
%
% INPUT
% showgraph             logical variable, set 1 to show plots of poincare sections
%
% EXAMPLE
% XXX

function [closedOrbitPosition, orbitPosition] = poincare_closed_orbit(flow,...
    vectorField,poincareSection,odeSolverOptions,showGraph)

narginchk(4,5)

if nargin == 4
    showGraph = true;
end

% Initial positions for Poincare orbits
orbitInitialPositionX = linspace(poincareSection.endPosition(1,1),...
    poincareSection.endPosition(2,1),poincareSection.numPoints);
orbitInitialPositionY = linspace(poincareSection.endPosition(1,2),...
    poincareSection.endPosition(2,2),poincareSection.numPoints);
orbitInitialPosition = transpose([orbitInitialPositionX; ...
    orbitInitialPositionY]);
initialDx = mean(sqrt(diff(orbitInitialPositionX).^2 + diff(orbitInitialPositionY).^2));

flowDomain = flow.domain;
flowResolution = flow.resolution;
orbitPosition = cell(poincareSection.numPoints,1);

% integrate orbits
parfor idx = 1:poincareSection.numPoints
    orbitPosition{idx} = integrate_line_closed(poincareSection.integrationLength,...
        orbitInitialPosition(idx,:),flowDomain,flowResolution,...
        vectorField,poincareSection.endPosition,odeSolverOptions);
end

% final position of orbits
orbitFinalPosition = cellfun(@(position)position(end,:),orbitPosition,...
    'UniformOutput',false);
orbitFinalPosition = cell2mat(orbitFinalPosition);

xLength = sqrt(diff(poincareSection.endPosition(:,1))^2 ...
    + diff(poincareSection.endPosition(:,2))^2);
if showGraph
    hfigure = figure;
    hAxes = axes;
    set(hAxes,'parent',hfigure);
    set(hAxes,'nextplot','add')
    set(hAxes,'box','on');
    set(hAxes,'xgrid','on');
    set(hAxes,'ygrid','on');
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
    plot(hAxes,abs(s(:,1)),zeros(length(abs(s(:,1))),1),'k--', 'linewidth', 1);       
    h0 = plot(hAxes,abs(s(:,1)),t(:,1)-s(:,1),'b.-','markersize',7);
    title('Closed orbit detection - Poincare section - return distance');
    legend(h0(1),'No closed orbits');
    xlabel(hAxes,'s'); ylabel(hAxes,'p(s) - s');
end

% find zero crossings of poincare return map (linear interpolation)
[~,closedOrbitInitialPosition] = crossing(t(:,1) - s(:,1),s(:,1));

if isempty(closedOrbitInitialPosition)
    
    closedOrbitInitialPosition(1,1:2) = [NaN NaN];
    closedOrbitPosition = [NaN NaN];
    
else
    
    if showGraph
        [nClosedOrbit ~] = size(closedOrbitInitialPosition);
        h1 = plot(hAxes,abs(closedOrbitInitialPosition),zeros(1,nClosedOrbit),'r.', 'markersize',7);
        legend([h1(1) 0], 'Zero crossings', 'No valid closed orbits');
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
    
    % FILTER:
    % (1) Discard closed orbits if images of neighbor points (P(x)) do not fall onto poincare section
    %***********************
    % PARAMETERS
    alphaThresh = 1e-4;
    distThresh = initialDx
    %***********************
    [nClosedOrbit ~] = size(closedOrbitInitialPosition);
    for i=1:nClosedOrbit
        % find 2 neighbor points of zero crossing
        [orbitInitialPositionSorted, ix] = sort(orbitInitialPosition(:,1));
        indx10 = max(find( closedOrbitInitialPosition(i,1) > orbitInitialPositionSorted ));
        indx20 = min(find( closedOrbitInitialPosition(i,1) < orbitInitialPositionSorted ));
        indx1 = min( ix(indx10), ix(indx20));
        indx2 = max( ix(indx10), ix(indx20));
        if indx2 <= indx1 || abs(indx1-indx2) ~=1
            error('Selection of neighbor orbits failed.')
        end
        % check if image of neighbor points falls on poincare section
        % i.e. angle of vectors to endpoints of poincare section equals pi
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
        % (2)
        % neighbor points of zero crossing must have a small return distance
        if  any( abs(t([indx1 indx2],1)-s([indx1 indx2],1)) > distThresh)
            closedOrbitInitialPosition(i,:) = NaN;
        end
    end            
    % erase invalid closed orbits
    [iy ~]= find(isnan(closedOrbitInitialPosition));
    closedOrbitInitialPosition(unique(iy),:) = [];
    
    if ~isempty(closedOrbitInitialPosition)
        [nClosedOrbit ~] = size(closedOrbitInitialPosition);
        % integrate closed orbits
        parfor idx = 1:nClosedOrbit
            closedOrbitPosition{idx} = integrate_line_closed(poincareSection.integrationLength,...
                closedOrbitInitialPosition(idx,:),flowDomain,flowResolution,...
                vectorField,poincareSection.endPosition,odeSolverOptions);
        end
        
        % FILTER: select outermost closed orbit
        s1(:,1) = closedOrbitInitialPosition(:,1) - poincareSection.endPosition(1,1);
        s1(:,2) = closedOrbitInitialPosition(:,2) - poincareSection.endPosition(1,2);
        distR = sqrt(s1(:,1).^2 + s1(:,2).^2);
        % plot all valid closed orbits
        if showGraph
            h2 = plot(hAxes,distR,0,'r^', 'markersize', 15);
        end        
        % outermost = largest distance from 1st point of poincare section
        indR = find( distR == max(distR) );
        closedOrbitInitialPosition = closedOrbitInitialPosition(indR,:);
        closedOrbitPosition = closedOrbitPosition{indR}(:,:);
        
        if showGraph
            h3 = plot(hAxes,distR(indR),0,'o', 'color', [0 0.6 0],'markersize', 20);
            legend([h1(1) h2(1), h3(1)], 'Closed orbits', 'Valid closed orbits', 'Outermost valid closed orbit');
            drawnow
        end
        
    else
        closedOrbitInitialPosition(1,1:2) = [NaN NaN];
        closedOrbitPosition = [NaN NaN];
        if showGraph
            drawnow
        end
    end
end


function alpha = angle2vectors(v1, v2)
v1_norm = sqrt(v1(1)^2 + v1(2)^2);
v2_norm = sqrt(v2(1)^2 + v2(2)^2);
alpha = acos( (v1(1)*v2(1)+v1(2)*v2(2))/(v1_norm*v2_norm) );

