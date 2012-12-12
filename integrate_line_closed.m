function position = integrate_line_closed(timespan,...
    initialCondition,domain,flowResolution,vectorGrid,odeSolverOptions)
%INTEGRATE_LINE Integrate line in non orientable vector field.

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

odeSolverOptions = odeset(odeSolverOptions,...
    'outputFcn',@(t,position,flag)ode_output(t,position,flag,...
    previousVector,vectorInterpolant,domain,flowResolution,vectorGrid),...
    'events',@(t,position)ode_events(t,position,domain));
sol = ode113(@(time,position)odefun(time,position,domain,...
    flowResolution,vectorGrid,vectorInterpolant,previousVector),...
    timespan,transpose(initialCondition),odeSolverOptions);
position = transpose(sol.y);
% position = transpose(deval(sol,linspace(0,sol.x(end),100)));

% FIXME Integration with event detection should not produce NaN positions
% nor positions outside domain in the first place. Need to research event
% detection accuracy.
position = remove_nan(position);
position = remove_outside(position,domain);

function output = odefun(~,position,domain,flowResolution,vectorGrid,...
    vectorInterpolant,previousVector)

continuousInterpolant = is_element_with_orient_discont(position,...
    domain,flowResolution,vectorGrid);

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

function status = ode_output(~,position,flag,vector,vectorInterpolant,...
    domain,flowResolution,vectorGrid)

if nargin < 3 || isempty(flag)

    % Use the last position
    position = position(:,end);
    
    continuousInterpolant = is_element_with_orient_discont(position,...
        domain,flowResolution,vectorGrid);

    if ~isempty(continuousInterpolant)
        currentVector = [continuousInterpolant.x(position(2),position(1)) ...
            continuousInterpolant.y(position(2),position(1))];
    else
        currentVector = [vectorInterpolant.x(position(2),position(1)) ...
            vectorInterpolant.y(position(2),position(1))];
    end
    
    vector.value = sign(currentVector*transpose(vector.value))*currentVector;

else
    switch(flag)
        case 'init'
            
            continuousInterpolant = is_element_with_orient_discont(...
                position,domain,flowResolution,vectorGrid);
            if ~isempty(continuousInterpolant)
                vector.value = [continuousInterpolant.x(position(2),...
                    position(1)) continuousInterpolant.y(position(2),...
                    position(1))];
            else
                vector.value = [vectorInterpolant.x(position(2),position(1)) ...
                    vectorInterpolant.y(position(2),position(1))];
            end
            
    end
end

status = 0;

function [distance,isTerminal,direction] = ode_events(time,position,domain)

if time < .1
    isTerminal = false;
else
    isTerminal = true;
end

direction = 1;

if any(isnan(position))
    distance = 0;
    return
end

% distance = drectangle(position,domain(1,1),domain(1,2),domain(2,1),...
%     domain(2,2));
% distance = position(1) - 4.15;
% distance = distPointToLineSegment([4.15; -.25],[4.15; -.6],position);
q1 = [4.15; -.6];
q2 = [4.15; -.3];
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/164048
distance = det([q2-q1,position-q1])/norm(q2-q1);
% disp(['event function distance: ',num2str(distance)])

function continuousInterpolant = ...
    is_element_with_orient_discont(position,domain,resolution,vector)
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
if drectangle(position,domain(1,1),domain(1,2),domain(2,1),...
    domain(2,2)) > 0
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

position1 = [(idxX-1)*deltaX+xMin (idxY-1)*deltaY+yMin];
vector1 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];

% Corner 2: upper-left
idxX = idxX - 1;
% position2 = [(idxX-1)*deltaX (idxY-1)*deltaY];
vector2 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];
if vector1*transpose(vector2) < 0
    isDiscontinuous = true;
    vector2 = -vector2;
end

% Corner 3: lower-left
idxY = idxY - 1;
position3 = [(idxX-1)*deltaX+xMin (idxY-1)*deltaY+yMin];
vector3 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];
if vector1*transpose(vector3) < 0
    isDiscontinuous = true;
    vector3 = -vector3;
end

% Corner 4: lower-right
idxX = idxX + 1;
% position4 = [(idxX-1)*deltaX (idxY-1)*deltaY];
vector4 = [vectorX(idxY,idxX) vectorY(idxY,idxX)];
if vector1*transpose(vector4) < 0
    isDiscontinuous = true;
    vector4 = -vector4;
end

if isDiscontinuous
    positionX = [position3(1) position1(1)];
    positionY = [position3(2) position1(2)];
    
    continuousInterpolant.x = griddedInterpolant({positionY,positionX},...
        [vector3(1) vector4(1); vector2(1) vector1(1)]);
    continuousInterpolant.y = griddedInterpolant({positionY,positionX},...
        [vector3(2) vector4(2); vector2(2) vector1(2)]);
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

function position = remove_outside(position,domain)
% Remove all positions from the first position outside domain onward

xMinIdx = find(position(:,1) < domain(1,1),1);
if isempty(xMinIdx)
    xMinIdx = size(position,1) + 1;
end
xMaxIdx = find(position(:,1) > domain(1,2),1);
if isempty(xMaxIdx)
    xMaxIdx = size(position,1) + 1;
end
yMinIdx = find(position(:,2) < domain(2,1),1);
if isempty(yMinIdx)
    yMinIdx = size(position,1) + 1;
end
yMaxIdx = find(position(:,2) > domain(2,2),1);
if isempty(yMaxIdx)
    yMaxIdx = size(position,1) + 1;
end
outsideIdx = min([xMinIdx xMaxIdx yMinIdx yMaxIdx]);
position = position(1:outsideIdx-1,:);
