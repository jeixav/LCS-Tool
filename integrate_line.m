%integrate_line Integrate line in non orientable vector field.
%
% SYNTAX
% position = integrate_line(timespan,initialCondition,domain,flowResolution,flowPeriodicBc,vectorGrid,odeSolverOptions)
% position = integrate_line(timespan,initialCondition,domain,flowResolution,flowPeriodicBc,vectorGrid,odeSolverOptions,poincareSection)
% position = integrate_line(...,'checkDiscontinuity',checkDiscontinuity)
% [position,iEvent] = integrate_line(...)
%
% OUTPUT ARGUMENTS
% iEvent: 1 if poincareSection specified and integration stopped at the
% Poincare section (i.e. position is a closed orbit.) 2 if poincareSection
% specified and integration stopped at domain boundaries. 3 if integration
% stopped at NaN. 0 otherwise (i.e. event detection did not terminate
% integration.)

function [position,varargout] = integrate_line(timespan,initialCondition,domain,flowResolution,flowPeriodicBc,vectorGrid,odeSolverOptions,varargin)

p = inputParser;
addOptional(p,'poincareSection',[])
addParameter(p,'checkDiscontinuity',true,@(input)validateattributes(input,{'logical'},{'scalar'}))
parse(p,varargin{:})

poincareSection = p.Results.poincareSection;
checkDiscontinuity = p.Results.checkDiscontinuity;

tmp = initialize_ic_grid(flowResolution,domain);
tmp = reshape(tmp(:,1),fliplr(flowResolution));
positionX = tmp(1,:);

tmp = initialize_ic_grid(flowResolution,domain);
tmp = reshape(tmp(:,2),fliplr(flowResolution));
positionY = tmp(:,1);

vectorXGrid = reshape(vectorGrid(:,1),fliplr(flowResolution));
vectorYGrid = reshape(vectorGrid(:,2),fliplr(flowResolution));

vectorInterpolant.x = griddedInterpolant({positionY,positionX},vectorXGrid);
vectorInterpolant.y = griddedInterpolant({positionY,positionX},vectorYGrid);

if isnan(vectorInterpolant.x(initialCondition(2),initialCondition(1)))
    warning([mfilename,':initialConditionNaN'],'vectorGrid at initialCondition = (%g,%g) is NaN',initialCondition(1),initialCondition(2))
    position = [NaN,NaN];
    return
end

previousVector = valueHandle;
previousVector.value = [];

discontinuousLargeAngle = valueHandle;
discontinuousLargeAngle.value = [];

if ~isempty(poincareSection)
    poincareSection = varargin{1};
    
    % set direction for event detection
    q1 = poincareSection(1,:);
    q2 = poincareSection(2,:);
    
    % vector field at initial position
    if checkDiscontinuity
        continuousInterpolant = is_element_with_orient_discont(initialCondition,domain,flowResolution,vectorGrid);
    end
    
    if checkDiscontinuity && ~isempty(continuousInterpolant)
        v(1) = continuousInterpolant.x([initialCondition(2),initialCondition(1)]);
        v(2) = continuousInterpolant.y([initialCondition(2),initialCondition(1)]);
    else
        v(1) = vectorInterpolant.x(initialCondition(2),initialCondition(1));
        v(2) = vectorInterpolant.y(initialCondition(2),initialCondition(1));
    end
    % vector along poincare section
    vPS = q2' - q1';
    dir = cross([vPS;0],[v';0]);
    if dir(3) > 0
        % look for zero crossing on rising edge
        direction = 1;
    else
        % look for zero crossing on falling edge
        direction = -1;
    end
end

odeSolverOptions = odeset(odeSolverOptions,'outputFcn',@(time,position,flag)ode_output(time,position,flag,previousVector,vectorInterpolant,domain,flowResolution,flowPeriodicBc,vectorGrid,checkDiscontinuity));
if isempty(poincareSection)
    odeSolverOptions = odeset(odeSolverOptions,'events',@(time,position)ode_events(time,position,domain,flowPeriodicBc));
else
    % FIXME Have not established poincareTimeTol = .1 is always suitable
    poincareTimeTol = .1;
    odeSolverOptions = odeset(odeSolverOptions,'events',@(time,position)ode_events_poincare(time,position,poincareSection,direction,domain,flowPeriodicBc,poincareTimeTol));
end

[~,position,tEvent,~,iEvent] = ode45(@(time,position)odefun(time,position,domain,flowResolution,flowPeriodicBc,vectorGrid,vectorInterpolant,previousVector,checkDiscontinuity,discontinuousLargeAngle),timespan,transpose(initialCondition),odeSolverOptions);

if nargout >= 2

    if ~isempty(poincareSection)
        % Discard event detections that occur before poincareTimeTol
        iEvent = iEvent(tEvent >= poincareTimeTol);
    end
    
    if isempty(iEvent)
        warning([mfilename,':iEventEmpty'],['iEvent empty; timespan may have to be increased. timespan = ',num2str(timespan(2))])
        iEvent = 0;
    end
    
    varargout{1} = iEvent(end);
end

if nargout == 3
    if ~isempty(discontinuousLargeAngle.value)
        varargout{2} = discontinuousLargeAngle.value;
    else
        varargout{2} = [];
    end
end

% FIXME Integration with event detection should not produce NaN positions
position = remove_nan(position);

function output = odefun(~,position,domain,flowResolution,flowPeriodicBc,vectorGrid,vectorInterpolant,previousVector,checkDiscontinuity,discontinuousLargeAngle)

position = transpose(apply_periodic_bc(transpose(position),flowPeriodicBc,domain));

if checkDiscontinuity
    [continuousInterpolant,lDiscontinuousLargeAngle] = is_element_with_orient_discont(position,domain,flowResolution,vectorGrid);
    % Record angle closest to pi/2 only
    if ~isempty(lDiscontinuousLargeAngle)
        if ~isempty(discontinuousLargeAngle.value)
            if abs(lDiscontinuousLargeAngle - pi/2) < abs(discontinuousLargeAngle.value - pi/2)
                discontinuousLargeAngle.value = lDiscontinuousLargeAngle;
            end
        else
            discontinuousLargeAngle.value = lDiscontinuousLargeAngle;
        end
    end
end

% ODE integrators expect column arrays
position = transpose(position);

if checkDiscontinuity && ~isempty(continuousInterpolant)
    output(:,1) = continuousInterpolant.x(position(:,2),position(:,1));
    output(:,2) = continuousInterpolant.y(position(:,2),position(:,1));
else
    output(:,1) = vectorInterpolant.x(position(:,2),position(:,1));
    output(:,2) = vectorInterpolant.y(position(:,2),position(:,1));
end

output = transpose(output);

% Orientation discontinuity
if ~isempty(previousVector.value) && ~all(isnan(previousVector.value))
    output = sign(previousVector.value*output)*output;
end

function status = ode_output(~,position,flag,vector,vectorInterpolant,domain,flowResolution,flowPeriodicBc,vectorGrid,checkDiscontinuity)

if nargin < 3 || isempty(flag)
    % Use last position
    position = position(:,end);
    position = transpose(apply_periodic_bc(transpose(position),flowPeriodicBc,domain));
    
    if checkDiscontinuity
        continuousInterpolant = is_element_with_orient_discont(position,domain,flowResolution,vectorGrid);
    end

    if checkDiscontinuity && ~isempty(continuousInterpolant)
        currentVector = [continuousInterpolant.x(position(2),position(1)),continuousInterpolant.y(position(2),position(1))];
    else
        currentVector = [vectorInterpolant.x(position(2),position(1)),vectorInterpolant.y(position(2),position(1))];
    end
    
    vector.value = sign(currentVector*transpose(vector.value))*currentVector;
else
    switch(flag)
        case 'init'
            if checkDiscontinuity
                continuousInterpolant = is_element_with_orient_discont(position,domain,flowResolution,vectorGrid);
            end
            
            if checkDiscontinuity && ~isempty(continuousInterpolant)
                vector.value = [continuousInterpolant.x(position(2),position(1)),continuousInterpolant.y(position(2),position(1))];
            else
                vector.value = [vectorInterpolant.x(position(2),position(1)),vectorInterpolant.y(position(2),position(1))];
            end
    end
end

status = 0;

function [distance,isTerminal,direction] = ode_events(~,position,domain,flowPeriodicBc)

isTerminal = true;
direction = 1;

if any(isnan(position))
    distance = 0;
    return
end
% shortest distance of position to domain boundaries
distance = drectangle(position,domain(1,1),domain(1,2),domain(2,1),domain(2,2),flowPeriodicBc);

function [distance,isTerminal,direction] = ode_events_poincare(time,position,poincareSection,directionPoincare,domain,flowPeriodicBc,poincareTimeTol)

% End-points of poincare section
q1 = poincareSection(1,:);
q2 = poincareSection(2,:);

% http://www.mathworks.com/matlabcentral/newsreader/view_thread/164048
% cross product of vector q1--q2 and vector q1--position
% positive on one side of poincare section, negative on other side
distancePoincare = det([q2 - q1;transpose(position) - q1])/norm(q2 - q1);

% Check that position is inside domain
[distanceRectangle,isTerminalRectangle,directionRectangle] = ode_events(time,position,domain,flowPeriodicBc);

if any(isnan(position))
    distance = [0;distanceRectangle];
    isTerminal = [true;isTerminalRectangle];
    direction = [directionPoincare;directionRectangle];
    return
end

if time < poincareTimeTol
    isTerminalPoincare = false;    
else
    isTerminalPoincare = true;
end

distance = [distancePoincare;distanceRectangle];
isTerminal = [isTerminalPoincare;isTerminalRectangle];
direction = [directionPoincare;directionRectangle];

function [continuousInterpolant,discontinuousLargeAngle] = is_element_with_orient_discont(position,domain,resolution,vector)
% Determine if position is between grid points with an orientation
% discontinuity in the vector field. If yes, return interpolant with 
% discontinuity removed.

continuousInterpolant = [];
discontinuousLargeAngle = [];

% FIXME This function should not get called in the first place if position
% is NaN.
if all(isnan(position))
    return
end

% FIXME Check if this function ever gets called if position is outside 
% domain.
if drectangle(position,domain(1,1),domain(1,2),domain(2,1),domain(2,2),[false,false]) >= 0
    return
end

vectorX = reshape(vector(:,1),fliplr(resolution));
vectorY = reshape(vector(:,2),fliplr(resolution));

isDiscontinuous = false;

% Index of Corner 1: upper-right. Regular grid spacing assumed.
deltaX = diff(domain(1,:))/(double(resolution(1)) - 1);
xMin = domain(1,1);
idxX = ceil((position(1) - xMin)/deltaX) + 1;

deltaY = diff(domain(2,:))/(double(resolution(2)) - 1);
yMin = domain(2,1);
idxY = ceil((position(2) - yMin)/deltaY) + 1;

position1 = [(idxX-1)*deltaX+xMin,(idxY-1)*deltaY+yMin];
vector1 = [vectorX(idxY,idxX),vectorY(idxY,idxX)];

% Corner 2: upper-left
idxX = idxX - 1;
vector2 = [vectorX(idxY,idxX),vectorY(idxY,idxX)];
% http://www.mathworks.com/matlabcentral/newsreader/view_original/381952
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/151925
angle = atan2(norm(cross([vector1,0],[vector2,0])),dot(vector1,vector2));

% FIXME smallAngle should not be hard-coded here
smallAngle = deg2rad(45);

if (angle > pi/2 - smallAngle) && (angle < pi/2 + smallAngle)
    discontinuousLargeAngle = angle;
end
if angle > pi/2
    isDiscontinuous = true;
    vector2 = -vector2;
end

% Corner 3: lower-left
idxY = idxY - 1;
position3 = [(idxX-1)*deltaX+xMin,(idxY-1)*deltaY+yMin];
vector3 = [vectorX(idxY,idxX),vectorY(idxY,idxX)];
angle = atan2(norm(cross([vector1,0],[vector3,0])),dot(vector1,vector3));
if (angle > pi/2 - smallAngle) && (angle < pi/2 + smallAngle)
    if ~isempty(discontinuousLargeAngle) && (angle > discontinuousLargeAngle)
        discontinuousLargeAngle = angle;
    else
        discontinuousLargeAngle = angle;
    end
end
if angle > pi/2
    isDiscontinuous = true;
    vector3 = -vector3;
end

% Corner 4: lower-right
idxX = idxX + 1;
vector4 = [vectorX(idxY,idxX),vectorY(idxY,idxX)];
angle = atan2(norm(cross([vector1,0],[vector4,0])),dot(vector1,vector4));
if (angle > pi/2 - smallAngle) && (angle < pi/2 + smallAngle)
    if ~isempty(discontinuousLargeAngle) && (angle > discontinuousLargeAngle)
        discontinuousLargeAngle = angle;
    else
        discontinuousLargeAngle = angle;
    end
end
if angle > pi/2
    isDiscontinuous = true;
    vector4 = -vector4;
end

if isDiscontinuous
    positionX = [position3(1),position1(1)];
    positionY = [position3(2),position1(2)];
    
    continuousInterpolant.x = griddedInterpolant({positionY,positionX},[vector3(1),vector4(1);vector2(1),vector1(1)]);
    continuousInterpolant.y = griddedInterpolant({positionY,positionX},[vector3(2),vector4(2);vector2(2),vector1(2)]);
end

function position = remove_nan(position)
% Remove all positions from the first NaN position onward

xNanIdx = find(isnan(position(:,1)),1);
if isempty(xNanIdx)
    xNanIdx = size(position,1) + 1;
end
yNanIdx = find(isnan(position(:,2)),1);
if isempty(yNanIdx)
    yNanIdx = size(position,1) + 1;
end
nanIdx = min([xNanIdx yNanIdx]);
position = position(1:nanIdx-1,:);

function d = drectangle(p,x1,x2,y1,y2,periodicBc)
% Compute signed distance function for rectangle with corners (x1,y1), 
% (x2,y1), (x1,y2), (x2,y2).
%
% Copied from DistMesh, http://persson.berkeley.edu/distmesh/.
% Copyright (C) 2004-2012 Per-Olof Persson.

d = -min(min(min(-y1+p(2),y2-p(2)),-x1+p(1)),x2-p(1));

if periodicBc(1)
    d = -min(-y1+p(2),y2-p(2));
end

if periodicBc(2)
    warning('Periodic BC in y not programmed')
end
