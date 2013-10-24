%integrate_line Integrate line in non orientable vector field.
%
% SYNTAX
% position = integrate_line(timespan,initialCondition,domain,flowResolution,flowPeriodicBc,vectorGrid,odeSolverOptions)
% position = integrate_line(timespan,initialCondition,domain,flowResolution,flowPeriodicBc,vectorGrid,odeSolverOptions,poincareSection)

function position = integrate_line(timespan,initialCondition,domain,flowResolution,flowPeriodicBc,vectorGrid,odeSolverOptions,varargin)

narginchk(7,8)

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

previousVector = valueHandle;
previousVector.value = [];

if nargin == 8
    poincareSection = varargin{1};
    
    % set direction for event detection
    q1 = poincareSection(1,:);
    q2 = poincareSection(2,:);
    
    % vector field at initial position
    continuousInterpolant = is_element_with_orient_discont(initialCondition,domain,flowResolution,vectorGrid);
    if ~isempty(continuousInterpolant)
        v(1) = continuousInterpolant.x([initialCondition(2),initialCondition(1)]);
        v(2) = continuousInterpolant.y([initialCondition(2),initialCondition(1)]);
    else
        v(1) = vectorInterpolant.x(initialCondition(2),initialCondition(1));
        v(2) = vectorInterpolant.y(initialCondition(2),initialCondition(1));
    end
    % vector along poincare section
    vPS = q2' - q1';
    dir = cross([vPS;0], [v';0]);
    if dir(3) > 0
        % look for zero crossing on rising edge
        direction = 1;
    else
        % look for zero crossing on falling edge
        direction = -1;
    end
end

odeSolverOptions = odeset(odeSolverOptions,'outputFcn',@(time,position,flag)ode_output(time,position,flag,previousVector,vectorInterpolant,domain,flowResolution,flowPeriodicBc,vectorGrid));
if nargin == 7
    odeSolverOptions = odeset(odeSolverOptions,'events',@(time,position)ode_events(time,position,domain,flowPeriodicBc));
else
    odeSolverOptions = odeset(odeSolverOptions,'events',@(time,position)ode_events_poincare(time,position,poincareSection,direction,domain,flowPeriodicBc));
    odeSolverOptions = odeset(odeSolverOptions,'initialStep',1e-9);
end

[~,position] = ode45(@(time,position)odefun(time,position,domain,flowResolution,flowPeriodicBc,vectorGrid,vectorInterpolant,previousVector),timespan,transpose(initialCondition),odeSolverOptions);

% FIXME Integration with event detection should not produce NaN positions
position = remove_nan(position);

function output = odefun(~,position,domain,flowResolution,flowPeriodicBc,vectorGrid,vectorInterpolant,previousVector)

position = transpose(apply_periodic_bc(transpose(position),flowPeriodicBc,domain));

continuousInterpolant = is_element_with_orient_discont(position,domain,flowResolution,vectorGrid);

% ODE integrators expect column arrays
position = transpose(position);

if ~isempty(continuousInterpolant)
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

function status = ode_output(~,position,flag,vector,vectorInterpolant,domain,flowResolution,flowPeriodicBc,vectorGrid)

if nargin < 3 || isempty(flag)
    % Use last position
    position = position(:,end);
    if size(position,2) > 1
        disp(size(position,2) > 1)
    end
    position = transpose(apply_periodic_bc(transpose(position),flowPeriodicBc,domain));
    
    continuousInterpolant = is_element_with_orient_discont(position,domain,flowResolution,vectorGrid);

    if ~isempty(continuousInterpolant)
        currentVector = [continuousInterpolant.x(position(2),position(1)),continuousInterpolant.y(position(2),position(1))];
    else
        currentVector = [vectorInterpolant.x(position(2),position(1)),vectorInterpolant.y(position(2),position(1))];
    end
    
    vector.value = sign(currentVector*transpose(vector.value))*currentVector;
else
    switch(flag)
        case 'init'
            continuousInterpolant = is_element_with_orient_discont(position,domain,flowResolution,vectorGrid);
            if ~isempty(continuousInterpolant)
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

function [distance,isTerminal,direction] = ode_events_poincare(time,position,poincareSection,direction,domain,flowPeriodicBc)

% end points of poincare section
q1 = poincareSection(1,:);
q2 = poincareSection(2,:);

% http://www.mathworks.com/matlabcentral/newsreader/view_thread/164048
% cross product of vector q1--q2 and vector q1--position
% positive on one side of poincare section, negative on other side
distancePoincare = det([q2 - q1;transpose(position) - q1])/norm(q2 - q1);

% Check distance to domain boundaries and then detect whether have crossed
% Poincare section or domain boundary
distanceRectangle = drectangle(position,domain(1,1),domain(1,2),domain(2,1),domain(2,2),flowPeriodicBc);

if abs(distancePoincare) < abs(distanceRectangle)
    distance = distancePoincare;
else
    distance = distanceRectangle;
    direction = 1;
end 

% isterminal = 1 if the integration is to terminate at a zero of this event function, otherwise, 0.
% FIXME Have not established time < .1 works in all cases.
if time < .1
    isTerminal = false;    
else
    isTerminal = true;
end

if any(isnan(position))
    distance = 0;
    return
end

function continuousInterpolant = is_element_with_orient_discont(position,domain,resolution,vector)
% Determine if position is between grid points with an orientation
% discontinuity in the vector field. If yes, return interpolant with 
% discontinuity removed.

continuousInterpolant = [];

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
% position2 = [(idxX-1)*deltaX (idxY-1)*deltaY];
vector2 = [vectorX(idxY,idxX),vectorY(idxY,idxX)];
if vector1*transpose(vector2) < 0
    isDiscontinuous = true;
    vector2 = -vector2;
end

% Corner 3: lower-left
idxY = idxY - 1;
position3 = [(idxX-1)*deltaX+xMin,(idxY-1)*deltaY+yMin];
vector3 = [vectorX(idxY,idxX),vectorY(idxY,idxX)];
if vector1*transpose(vector3) < 0
    isDiscontinuous = true;
    vector3 = -vector3;
end

% Corner 4: lower-right
idxX = idxX + 1;
% position4 = [(idxX-1)*deltaX (idxY-1)*deltaY];
vector4 = [vectorX(idxY,idxX),vectorY(idxY,idxX)];
if vector1*transpose(vector4) < 0
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
