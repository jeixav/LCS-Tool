function noStretchLine = compute_no_stretch_line(flow,noStretchLine,verbose)

narginchk(2,3)

if nargin < 3
    verbose = false;
end

if ~all(isfield(noStretchLine,{'chiPos','chiNeg'}))
    % Calculate etaPos and etaNeg on the grid on which the flow is
    % calculated
    l1 = flow.cgEigenvalue(:,1);
    l2 = flow.cgEigenvalue(:,2);
    xi1 = flow.cgEigenvector(:,1:2);
    xi2 = flow.cgEigenvector(:,3:4);

    alpha = sqrt((l2-1)./(l2-l1));
    if ~isreal(alpha)
        warning('compute_no_stretch_line:alpha_not_real','alpha not real')
    end
    
    beta = sqrt((1-l1)./(l2-l1));
    if ~isreal(beta)
        warning('compute_no_stretch_line:beta_not_real','beta not real')
    end
    
    if ~isfield(noStretchLine,'chiPos')
        noStretchLine.chiPos(:,1) = alpha.*xi1(:,1) + beta.*xi2(:,1);
        noStretchLine.chiPos(:,2) = alpha.*xi1(:,2) + beta.*xi2(:,2);
    end
    
    if ~isfield(noStretchLine,'chiNeg')
        noStretchLine.chiNeg(:,1) = alpha.*xi1(:,1) - beta.*xi2(:,1);
        noStretchLine.chiNeg(:,2) = alpha.*xi1(:,2) - beta.*xi2(:,2);
    end
end

if isfield(noStretchLine,'resolution')
    noStretchLine.initialPosition = initialize_ic_grid(...
        noStretchLine.resolution,flow.domain);
end

if verbose
    progressBar = ConsoleProgressBar;
    progressBar.setText(mfilename)
    progressBar.setTextPosition('left')
    progressBar.setElapsedTimeVisible(1);
    progressBar.setRemainedTimeVisible(1);
    progressBar.setLength(20)
    progressBar.setMaximum(numel(noStretchLine.initialPosition))
    progressBar.start
    % FIXME ConsoleProgressBar does not work with SPMD
    % noStretchLine.odeSolverOptions = odeset(noStretchLine.odeSolverOptions,...
    %   'outputFcn',@(t,y,flag)progress_bar(t,y,flag,progressBar));
else
    progressBar = false;
end

if ~isfield(noStretchLine,'odeSolverOptions')
    noStretchLine.odeSolverOptions = [];
end

nNoStretchLines = size(noStretchLine.initialPosition,1);
iNoStretchLine = 1:nNoStretchLines;

timespan = [0 noStretchLine.finalTime];

spmdFun = @(idx)integrate_individual_line(timespan,...
    noStretchLine.initialPosition(idx,:),flow.domain,flow.resolution,...
    noStretchLine.chiPos,noStretchLine.odeSolverOptions);

spmd
    iNoStretchLineC = codistributed(iNoStretchLine);
    iNoStretchLineL = getLocalPart(iNoStretchLineC);
    positionPosL = arrayfun(spmdFun,iNoStretchLineL,'UniformOutput',false);
    positionPosC = codistributed.build(positionPosL,...
        getCodistributor(iNoStretchLineC));
end

noStretchLine.positionPos = gather(positionPosC);

spmdFun = @(idx)integrate_individual_line(timespan,...
    noStretchLine.initialPosition(idx,:),flow.domain,flow.resolution,...
    noStretchLine.chiNeg,noStretchLine.odeSolverOptions);

spmd
    iNoStretchLineC = codistributed(iNoStretchLine);
    iNoStretchLineL = getLocalPart(iNoStretchLineC);
    positionNegL = arrayfun(spmdFun,iNoStretchLineL,'UniformOutput',false);
    positionNegC = codistributed.build(positionNegL,...
        getCodistributor(iNoStretchLineC));
end

noStretchLine.positionNeg = gather(positionNegC);

% Backward time integration
timespan = -timespan;

spmdFun = @(idx)integrate_individual_line(timespan,...
    noStretchLine.initialPosition(idx,:),flow.domain,flow.resolution,...
    noStretchLine.chiPos,noStretchLine.odeSolverOptions);

spmd
    iNoStretchLineC = codistributed(iNoStretchLine);
    iNoStretchLineL = getLocalPart(iNoStretchLineC);
    positionPosL = arrayfun(spmdFun,iNoStretchLineL,'UniformOutput',false);
    positionPosC = codistributed.build(positionPosL,...
        getCodistributor(iNoStretchLineC));
end

spmdFun = @(idx)integrate_individual_line(timespan,...
    noStretchLine.initialPosition(idx,:),flow.domain,flow.resolution,...
    noStretchLine.chiNeg,noStretchLine.odeSolverOptions);

spmd
    iNoStretchLineC = codistributed(iNoStretchLine);
    iNoStretchLineL = getLocalPart(iNoStretchLineC);
    positionNegL = arrayfun(spmdFun,iNoStretchLineL,'UniformOutput',false);
    positionNegC = codistributed.build(positionNegL,...
        getCodistributor(iNoStretchLineC));
end

% Concatenate forward and backward time integration results
tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionPosC),...
    'uniformOutput',false);
noStretchLine.positionPos = cellfun(@(a,b)[a;b],tmp,...
    noStretchLine.positionPos,'UniformOutput',false);

tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionNegC),...
    'uniformOutput',false);
noStretchLine.positionNeg = cellfun(@(a,b)[a;b],tmp,...
    noStretchLine.positionNeg,'UniformOutput',false);

if verbose
    progressBar.setValue(progressBar.maximum)
    progressBar.stop
    fprintf('\n')
end

is_position_outside_domain(noStretchLine.positionPos,flow.domain)
is_position_outside_domain(noStretchLine.positionNeg,flow.domain)

function position = integrate_individual_line(timespan,...
    initialCondition,domain,flowResolution,vectorGrid,odeSolverOptions)

positionX = linspace(domain(1,1),domain(1,2),flowResolution(1));
positionY = linspace(domain(2,1),domain(2,2),flowResolution(2));
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
sol = ode45(@(time,position)odefun(time,position,domain,...
    flowResolution,vectorGrid,vectorInterpolant,previousVector),timespan,...
    transpose(initialCondition),odeSolverOptions);
position = transpose(sol.y);

% Remove NaN values of final integration step.
position = position(~isnan(position(:,1)),:);

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

function [distance,isTerminal,direction] = ode_events(~,position,domain)

isTerminal = true;
direction = 1;

if any(isnan(position))
    distance = 1;
    return
end

distance = drectangle(position,domain(1,1),domain(1,2),domain(2,1),...
    domain(2,2));

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
