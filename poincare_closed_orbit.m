% poincare_closed_orbit Find closed orbits using Poincare section return
% map
%
% SYNTAX
% [closedOrbitPosition,orbitPosition] = poincare_closed_orbit(flow,vectorField,poincareSection,odeSolverOptions,nBisection,dThresh,showGraph)
%
% INPUT ARGUMENTS
% showgraph: logical variable, set true to show plots of Poincare sections

function [closedOrbitPosition,orbitPosition] = poincare_closed_orbit(flow,vectorField,poincareSection,odeSolverOptions,nBisection,dThresh,showGraph)

narginchk(6,7)

if nargin == 6
    showGraph = false;
end

p = inputParser;
p.KeepUnmatched = true;

addParamValue(p,'periodicBc',[false,false],@(i)validateattributes(i,{'logical'},{'size',[1,2]}));

parse(p,flow)
flow.periodicBc = p.Results.periodicBc;

% Poincare section vector
p = poincareSection.endPosition(2,:) - poincareSection.endPosition(1,:);

% Initial positions for Poincare orbits
orbitInitialPositionX = linspace(poincareSection.endPosition(1,1),poincareSection.endPosition(2,1),poincareSection.numPoints);
orbitInitialPositionY = linspace(poincareSection.endPosition(1,2),poincareSection.endPosition(2,2),poincareSection.numPoints);
orbitInitialPosition = transpose([orbitInitialPositionX;orbitInitialPositionY]);

orbitPosition = cell(poincareSection.numPoints,1);

% integrate orbits
parfor idx = 1:poincareSection.numPoints
    orbitPosition{idx} = integrate_line(poincareSection.integrationLength,orbitInitialPosition(idx,:),flow.domain,flow.resolution,flow.periodicBc,vectorField,odeSolverOptions,poincareSection.endPosition); %#ok<PFBNS>
end

% final position of orbits
orbitFinalPosition = cellfun(@(position)position(end,:),orbitPosition,'UniformOutput',false);
orbitFinalPosition = cell2mat(orbitFinalPosition);

xLength = hypot(diff(poincareSection.endPosition(:,1)),diff(poincareSection.endPosition(:,2)));

% angle of Poincare section with x-axis
theta = -atan((poincareSection.endPosition(1,2) - poincareSection.endPosition(2,2))/(poincareSection.endPosition(1,1) - poincareSection.endPosition(2,1)));
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
    hfigure = figure;
    hAxes = axes;
    set(hAxes,'parent',hfigure)
    set(hAxes,'nextplot','add')
    title(hAxes,'Poincare return map')
    set(hAxes,'box','on')
    set(hAxes,'xgrid','on')
    set(hAxes,'ygrid','on')
    set(hAxes,'xlim',[0,xLength])
    hReturnMap = plot(hAxes,abs(s(:,1)),t(:,1)-s(:,1));
    set(hReturnMap,'color','b')
    xlabel(hAxes,'s')
    ylabel(hAxes,'p(s) - s')
end

% find zero crossings of Poincare return map (linear interpolation)
[~,closedOrbitInitialPosition] = crossing(t(:,1) - s(:,1),s(:,1));

if isempty(closedOrbitInitialPosition)
    closedOrbitPosition{1} = [NaN,NaN];
else
    if showGraph
        hClosedOrbitInitialPosition = plot(hAxes,abs(closedOrbitInitialPosition),zeros(size(closedOrbitInitialPosition)));
        set(hClosedOrbitInitialPosition,'LineStyle','none')
        set(hClosedOrbitInitialPosition,'marker','o')
        set(hClosedOrbitInitialPosition,'MarkerEdgeColor','b')
        set(hClosedOrbitInitialPosition,'DisplayName','Zero crossing candidate')
        hLegend = legend(hAxes,hClosedOrbitInitialPosition);
        set(hLegend,'location','best')
        drawnow
    end
    
    % Rotate to theta
    xx = [transpose(closedOrbitInitialPosition),zeros(numel(closedOrbitInitialPosition),1)];
    xx = rotationMatrix\transpose(xx);
    xx = transpose(xx);
    
    % Translate from origin
    closedOrbitInitialPositionX = xx(:,1) + poincareSection.endPosition(1,1);
    closedOrbitInitialPositionY = xx(:,2) + poincareSection.endPosition(1,2);
    
    closedOrbitInitialPosition = [closedOrbitInitialPositionX,closedOrbitInitialPositionY];
    
    % FILTER:
    % Discard discontinuous zero crossings
    % Refine zero crossings with bisection method
    % PARAMETERS
    distThresh = dThresh * xLength;

    nZeroCrossing = size(closedOrbitInitialPosition,1);
    for i = 1:nZeroCrossing
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
        
        for j = 1:nBisection
            % get return distance for p1, p2
            p1finalPos = integrate_line(poincareSection.integrationLength,p1,flow.domain,flow.resolution,flow.periodicBc,vectorField,odeSolverOptions,poincareSection.endPosition);
            p1end = p1finalPos(end,:);
            p1dist = dot(p1end - p1,p/norm(p));
            p2finalPos = integrate_line(poincareSection.integrationLength,p2,flow.domain,flow.resolution,flow.periodicBc,vectorField,odeSolverOptions,poincareSection.endPosition);
            p2end = p2finalPos(end,:);
            p2dist = dot(p2end - p2,p/norm(p));
            
            % bisect
            p3 = (p1+p2)/2;            
            % return distance for p3
            p3finalPos = integrate_line(poincareSection.integrationLength,p3,flow.domain,flow.resolution,flow.periodicBc,vectorField,odeSolverOptions,poincareSection.endPosition);
            p3end = p3finalPos(end,:);
            p3dist = dot(p3end - p3,p/norm(p));
            
            if j ~= nBisection
                if p1dist*p3dist < 0
                    p2 = p3;
                else
                    p1 = p3;
                end
            end
        end
        % neighbor points of zero crossing must have a small return distance
        if any(abs([p1dist,p2dist]) > distThresh)
            closedOrbitInitialPosition(i,:) = NaN;
        else
            closedOrbitInitialPosition(i,:) = p3;
        end
    end
    % Erase invalid closed orbits
    [iy,~] = find(isnan(closedOrbitInitialPosition));
    closedOrbitInitialPosition(unique(iy),:) = [];
    
    if ~isempty(closedOrbitInitialPosition)
        nClosedOrbit = size(closedOrbitInitialPosition,1);
        % Integrate closed orbits
        parfor idx = 1:nClosedOrbit
            closedOrbitPosition{idx} = integrate_line(poincareSection.integrationLength,closedOrbitInitialPosition(idx,:),flow.domain,flow.resolution,flow.periodicBc,vectorField,odeSolverOptions,poincareSection.endPosition); %#ok<PFBNS>
        end
        
        % FILTER: select outermost closed orbit
        s1(:,1) = closedOrbitInitialPosition(:,1) - poincareSection.endPosition(1,1);
        s1(:,2) = closedOrbitInitialPosition(:,2) - poincareSection.endPosition(1,2);
        distR = hypot(s1(:,1),s1(:,2));

        % Plot all valid zero crossings
        if showGraph
            delete(hLegend)
            delete(hClosedOrbitInitialPosition)
            hZeroCrossing = plot(hAxes,distR,zeros(size(distR)));
            set(hZeroCrossing,'LineStyle','none')
            set(hZeroCrossing,'marker','o')
            set(hZeroCrossing,'MarkerEdgeColor','b')
            set(hZeroCrossing,'MarkerFaceColor','b')
            set(hZeroCrossing,'DisplayName','Zero crossing')
            hLegend = legend(hAxes,hZeroCrossing);
            set(hLegend,'location','best')
            drawnow
        end
        
        % Sort closed orbits: closedOrbitPosition{1} is innermost,
        % closedOrbitPosition{end} is outermost.
        [~,sortIndex] = sort(distR);
        closedOrbitPosition = closedOrbitPosition(sortIndex);
        
        % Plot outermost zero crossing
        if showGraph
            hOutermostZeroCrossing = plot(hAxes,distR(sortIndex(end)),0);
            set(hOutermostZeroCrossing,'LineStyle','none')
            set(hOutermostZeroCrossing,'marker','o')
            set(hOutermostZeroCrossing,'MarkerEdgeColor','r')
            set(hOutermostZeroCrossing,'MarkerFaceColor','r')
            set(hOutermostZeroCrossing,'DisplayName','Outermost zero crossing')
            hLegend = legend(hAxes,[hZeroCrossing;hOutermostZeroCrossing]);
            set(hLegend,'Location','best')
            drawnow
        end
        
    else
        closedOrbitPosition{1} = [NaN,NaN];
    end
end

function [ind,t0,s0] = crossing(S,t,level,imeth)
% CROSSING find the crossings of a given level of a signal
%   ind = CROSSING(S) returns an index vector ind, the signal
%   S crosses zero at ind or at between ind and ind+1
%   [ind,t0] = CROSSING(S,t) additionally returns a time
%   vector t0 of the zero crossings of the signal S. The crossing
%   times are linearly interpolated between the given times t
%   [ind,t0] = CROSSING(S,t,level) returns the crossings of the
%   given level instead of the zero crossings
%   ind = CROSSING(S,[],level) as above but without time interpolation
%   [ind,t0] = CROSSING(S,t,level,par) allows additional parameters
%   par = {'none'|'linear'}.
%	With interpolation turned off (par = 'none') this function always
%	returns the value left of the zero (the data point thats nearest
%   to the zero AND smaller than the zero crossing).
%
%	[ind,t0,s0] = ... also returns the data vector corresponding to 
%	the t0 values.
%
%	[ind,t0,s0,t0close,s0close] additionally returns the data points
%	closest to a zero crossing in the arrays t0close and s0close.
%
%	This version has been revised incorporating the good and valuable
%	bugfixes given by users on Matlabcentral. Special thanks to
%	Howard Fishman, Christian Rothleitner, Jonathan Kellogg, and
%	Zach Lewis for their input. 

% Steffen Brueckner, 2002-09-25
% Steffen Brueckner, 2007-08-27		revised version

% Copyright (c) Steffen Brueckner, 2002-2007
% brueckner@sbrs.net

% Source: http://www.mathworks.com/matlabcentral/fileexchange/2432

% check the number of input arguments
narginchk(1,4);

% check the time vector input for consistency
if nargin < 2 || isempty(t)
	% if no time vector is given, use the index vector as time
    t = 1:length(S);
elseif length(t) ~= length(S)
	% if S and t are not of the same length, throw an error
    error('t and S must be of identical length!');    
end

% check the level input
if nargin < 3
	% set standard value 0, if level is not given
    level = 0;
end

% check interpolation method input
if nargin < 4
    imeth = 'linear';
end

% make row vectors
t = t(:)';
S = S(:)';

% always search for zeros. So if we want the crossing of 
% any other threshold value "level", we subtract it from
% the values and search for zeros.
S   = S - level;

% first look for exact zeros
ind0 = find( S == 0 ); 

% then look for zero crossings between data points
S1 = S(1:end-1) .* S(2:end);
ind1 = find( S1 < 0 );

% bring exact zeros and "in-between" zeros together 
ind = sort([ind0 ind1]);

% and pick the associated time values
t0 = t(ind); 
s0 = S(ind);

if strcmp(imeth,'linear')
    % linear interpolation of crossing
    for ii=1:length(t0)
        if abs(S(ind(ii))) > eps(S(ind(ii)))
            % interpolate only when data point is not already zero
            NUM = (t(ind(ii)+1) - t(ind(ii)));
            DEN = (S(ind(ii)+1) - S(ind(ii)));
            DELTA =  NUM / DEN;
            t0(ii) = t0(ii) - S(ind(ii)) * DELTA;
            % I'm a bad person, so I simply set the value to zero
            % instead of calculating the perfect number ;)
            s0(ii) = 0;
        end
    end
end

% FIXME Causes error when minimum is between first and second point
% % Addition:
% % Some people like to get the data points closest to the zero crossing,
% % so we return these as well
% [~,II] = min(abs([S(ind-1) ; S(ind) ; S(ind+1)]),[],1); 
% ind2 = ind + (II-2); %update indices 
% 
% t0close = t(ind2);
% s0close = S(ind2);
