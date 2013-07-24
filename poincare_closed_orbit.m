% poincare_closed_orbit Find closed orbits using Poincare section return
% map
%
% SYNTAX
%[closedOrbitPosition,orbitPosition] = poincare_closed_orbit(flow,vectorField,poincareSection,odeSolverOptions,showGraph)
%
% INPUT ARGUMENTS
% showgraph: logical variable, set true to show plots of Poincare sections

function [closedOrbitPosition,orbitPosition] = poincare_closed_orbit(flow,...
    vectorField,poincareSection,odeSolverOptions,nBisection,dThresh,showGraph)

narginchk(6,7)

if nargin == 6
    showGraph = true;
end

% poincare section vector
p = poincareSection.endPosition(2,:) - poincareSection.endPosition(1,:);

% Initial positions for Poincare orbits
orbitInitialPositionX = linspace(poincareSection.endPosition(1,1),...
    poincareSection.endPosition(2,1),poincareSection.numPoints);
orbitInitialPositionY = linspace(poincareSection.endPosition(1,2),...
    poincareSection.endPosition(2,2),poincareSection.numPoints);
orbitInitialPosition = transpose([orbitInitialPositionX; ...
    orbitInitialPositionY]);

flowDomain = flow.domain;
flowResolution = flow.resolution;
orbitPosition = cell(poincareSection.numPoints,1);

% integrate orbits
parfor idx = 1:poincareSection.numPoints
    orbitPosition{idx} = integrate_line_closed(poincareSection.integrationLength,...
        orbitInitialPosition(idx,:),flowDomain,flowResolution,...
        vectorField,poincareSection.endPosition,odeSolverOptions); %#ok<PFBNS>
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
    
    closedOrbitInitialPosition(1,1:2) = [NaN NaN]; %#ok<NASGU>
    closedOrbitPosition = [NaN NaN];
    
else
    
    if showGraph
        nClosedOrbit = size(closedOrbitInitialPosition,1);
        h1 = plot(hAxes,abs(closedOrbitInitialPosition),zeros(1,nClosedOrbit),'r.', 'markersize',7);
        legend([h1(1) 0], 'Zero crossings', 'No valid closed orbits');
        drawnow
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
    % Discard discontinuous zero crossings
    % Refine zero crossings with bisection method
    %***********************
    % PARAMETERS
    distThresh = dThresh * xLength;
    %***********************
    nZeroCrossing = size(closedOrbitInitialPosition,1);
    for i=1:nZeroCrossing
        % find 2 neighbor points of zero crossing
        [orbitInitialPositionSorted, ix] = sort(orbitInitialPosition(:,1));
        indx10 = find(closedOrbitInitialPosition(i,1) > orbitInitialPositionSorted,1,'last');
        indx20 = find(closedOrbitInitialPosition(i,1) < orbitInitialPositionSorted,1,'first');
        indx1 = min( ix(indx10), ix(indx20));
        indx2 = max( ix(indx10), ix(indx20));
        if indx2 <= indx1 || abs(indx1-indx2) ~=1
            error('Selection of neighbor orbits failed.')
        end        
        % Bisection method
        % neighbor points
        p1 = orbitInitialPosition(indx1,:);
        p2 = orbitInitialPosition(indx2,:);
        
        for j=1:nBisection
            % get return distance for p1, p2
            p1finalPos = integrate_line_closed(poincareSection.integrationLength,...
                p1,flowDomain,flowResolution,vectorField,poincareSection.endPosition,...
                odeSolverOptions);
            p1end = p1finalPos(end,:);
            p1dist = dot(p1end - p1,p/norm(p));
            p2finalPos = integrate_line_closed(poincareSection.integrationLength,...
                p2,flowDomain,flowResolution,vectorField,poincareSection.endPosition,...
                odeSolverOptions);
            p2end = p2finalPos(end,:);
            p2dist = dot(p2end - p2,p/norm(p));
            
            % bisect
            p3 = (p1+p2)/2;            
            % return distance for p3
            p3finalPos = integrate_line_closed(poincareSection.integrationLength,...
                p3,flowDomain,flowResolution,...
                vectorField,poincareSection.endPosition,odeSolverOptions);
            p3end = p3finalPos(end,:);
            p3dist = dot(p3end - p3,p/norm(p));
            
            % plot - check bisection method
            if showGraph
                plot(hAxes, ...
                    [dot(p1-poincareSection.endPosition(1,:),p/norm(p)) ...
                    dot(p3-poincareSection.endPosition(1,:),p/norm(p)) ...
                    dot(p2-poincareSection.endPosition(1,:),p/norm(p))],...
                    [p1dist p3dist p2dist],'k-+');
            end            
            if j~=nBisection
                if p1dist*p3dist < 0
%                     p1 = p1;
                    p2 = p3;
                else
                    p1 = p3;
%                     p2 = p2;
                end
            end
        end
        % neighbor points of zero crossing must have a small return distance
        if  any( abs([p1dist p2dist]) > distThresh)
            closedOrbitInitialPosition(i,:) = NaN;
        else
            closedOrbitInitialPosition(i,:) = p3;
        end
    end
    % Erase invalid closed orbits
    [iy,~]= find(isnan(closedOrbitInitialPosition));
    closedOrbitInitialPosition(unique(iy),:) = [];
    
    if ~isempty(closedOrbitInitialPosition)
        nClosedOrbit = size(closedOrbitInitialPosition,1);
        % Integrate closed orbits
        parfor idx = 1:nClosedOrbit
            closedOrbitPosition{idx} = integrate_line_closed(poincareSection.integrationLength,...
                closedOrbitInitialPosition(idx,:),flowDomain,flowResolution,...
                vectorField,poincareSection.endPosition,odeSolverOptions); %#ok<PFBNS>
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
        closedOrbitInitialPosition = closedOrbitInitialPosition(indR,:); %#ok<NASGU>
        closedOrbitPosition = closedOrbitPosition{indR}(:,:);
        
        if showGraph
            h3 = plot(hAxes,distR(indR),0,'o', 'color', [0 0.6 0],'markersize', 20);
            legend([h1(1) h2(1), h3(1)], 'Closed orbits', 'Valid closed orbits', 'Outermost valid closed orbit');
            drawnow
        end
        
    else
        closedOrbitInitialPosition(1,1:2) = [NaN NaN]; %#ok<NASGU>
        closedOrbitPosition = [NaN NaN];
        if showGraph
            drawnow
        end
    end
end

