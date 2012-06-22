function shearline = compute_shearline(flow,shearline,verbose)

narginchk(2,3)

if nargin < 3
    verbose = false;
end

if ~all(isfield(shearline,{'etaPos','etaNeg'}))
    % Calculate etaPos and etaNeg on the grid on which the flow is
    % calculated
    l1 = flow.cgEigenvalue(:,1);
    l2 = flow.cgEigenvalue(:,2);
    xi1 = flow.cgEigenvector(:,1:2);
    xi2 = flow.cgEigenvector(:,3:4);
    
    if ~isfield(shearline,'etaPos')
        shearline.etaPos(:,1) = sqrt(sqrt(l2)./(sqrt(l1) ...
            + sqrt(l2))).*xi1(:,1) + sqrt(sqrt(l1)./(sqrt(l1) ...
            + sqrt(l2))).*xi2(:,1);
        shearline.etaPos(:,2) = sqrt(sqrt(l2)./(sqrt(l1) ...
            + sqrt(l2))).*xi1(:,2) + sqrt(sqrt(l1)./(sqrt(l1) ...
            + sqrt(l2))).*xi2(:,2);
        if ~isreal(shearline.etaPos)
            warning('compute_shearline:etaPos_not_real','etaPos not real')
        end
    end
    
    if ~isfield(shearline,'etaNeg')
        shearline.etaNeg(:,1) = sqrt(sqrt(l2)./(sqrt(l1) ...
            + sqrt(l2))).*xi1(:,1) - sqrt(sqrt(l1)./(sqrt(l1) ...
            + sqrt(l2))).*xi2(:,1);
        shearline.etaNeg(:,2) = sqrt(sqrt(l2)./(sqrt(l1) ...
            + sqrt(l2))).*xi1(:,2) - sqrt(sqrt(l1)./(sqrt(l1) ...
            + sqrt(l2))).*xi2(:,2);
        if ~isreal(shearline.etaNeg)
            warning('compute_shearline:etaNeg_not_real','etaNeg not real')
        end
    end
end

timespan = [0 shearline.finalTime];

if isfield(shearline,'resolution')
    shearline.initialPosition = initialize_ic_grid(shearline.resolution,...
        flow.domain);
end

if verbose
    progressBar = ConsoleProgressBar;
    progressBar.setText(mfilename)
    progressBar.setTextPosition('left')
    progressBar.setElapsedTimeVisible(1);
    progressBar.setRemainedTimeVisible(1);
    progressBar.setLength(20)
    progressBar.setMaximum(numel(shearline.initialPosition))
    progressBar.start
    % FIXME ConsoleProgressBar does not work with SPMD
    % shearline.odeSolverOptions = odeset(shearline.odeSolverOptions,...
    %   'outputFcn',@(t,y,flag)progress_bar(t,y,flag,progressBar));
else
    progressBar = false;
end

nShearlines = size(shearline.initialPosition,1);
iShearline = 1:nShearlines;

if ~isfield(shearline,'odeSolverOptions')
    shearline.odeSolverOptions = [];
end

spmdFunPos = @(idx)integrate_individual_shearline(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaPos,shearline.odeSolverOptions);
spmdFunNeg = @(idx)integrate_individual_shearline(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaNeg,shearline.odeSolverOptions);
spmd
    iShearlineC = codistributed(iShearline);
    iShearlineL = getLocalPart(iShearlineC);
    positionPosL = arrayfun(spmdFunPos,iShearlineL,'UniformOutput',false);
    positionNegL = arrayfun(spmdFunNeg,iShearlineL,'UniformOutput',false);
    positionPosC = codistributed.build(positionPosL,...
        getCodistributor(iShearlineC));
    positionNegC = codistributed.build(positionNegL,...
        getCodistributor(iShearlineC));
end
shearline.positionPos = gather(positionPosC);
shearline.positionNeg = gather(positionNegC);

% Backward time integration
timespan = -timespan;
spmdFunPos = @(idx)integrate_individual_shearline(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaPos,shearline.odeSolverOptions);
spmdFunNeg = @(idx)integrate_individual_shearline(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaNeg,shearline.odeSolverOptions);
spmd
    iShearlineC = codistributed(iShearline);
    iShearlineL = getLocalPart(iShearlineC);
    positionPosL = arrayfun(spmdFunPos,iShearlineL,'UniformOutput',false);
    positionNegL = arrayfun(spmdFunNeg,iShearlineL,'UniformOutput',false);
    positionPosC = codistributed.build(positionPosL,...
        getCodistributor(iShearlineC));
    positionNegC = codistributed.build(positionNegL,...
        getCodistributor(iShearlineC));
end

% Concatenate forward and backward time integration results
tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionPosC),...
    'uniformOutput',false);
shearline.positionPos = cellfun(@(a,b)[a;b],tmp,shearline.positionPos,...
    'UniformOutput',false);
tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionNegC),...
    'uniformOutput',false);
shearline.positionNeg = cellfun(@(a,b)[a;b],tmp,shearline.positionNeg,...
    'UniformOutput',false);

if verbose
    progressBar.setValue(progressBar.maximum)
    progressBar.stop
    fprintf('\n')
end

function position = integrate_individual_shearline(timespan,shearlineIc,...
    domain,flowResolution,etaGrid,odeSolverOptions)

positionX = linspace(domain(1,1),domain(1,2),flowResolution(1));
positionY = linspace(domain(2,1),domain(2,2),flowResolution(2));
etaXGrid = reshape(etaGrid(:,1),fliplr(flowResolution));
etaYGrid = reshape(etaGrid(:,2),fliplr(flowResolution));

etaInterpolant.X = griddedInterpolant({positionY,positionX},etaXGrid);
etaInterpolant.Y = griddedInterpolant({positionY,positionX},etaYGrid);

previousEta = valueHandle;
previousEta.value = [];

odeSolverOptions = odeset(odeSolverOptions,...
    'outputFcn',@(t,position,flag)ode_output(t,position,flag,...
    previousEta,etaInterpolant,domain,flowResolution,etaGrid),...
    'events',@(t,position)ode_events(t,position,domain));
sol = ode45(@(time,position)odefun(time,position,domain,...
    flowResolution,etaGrid,etaInterpolant,previousEta),timespan,...
    transpose(shearlineIc),odeSolverOptions);
position = transpose(sol.y);

% FIXME Integration with event detection should not produce NaN positions
% nor positions outside domain in the first place. Need to research event
% detection accuracy.
position = remove_nan(position);
position = remove_outside(position,domain);

function output = odefun(~,position,domain,flowResolution,etaGrid,...
    etaInterpolant,previousEta)

continuousInterpolant = is_element_with_orient_discont(position,...
    domain,flowResolution,etaGrid);

% ODE integrators expect column arrays
position = transpose(position);

if ~isempty(continuousInterpolant)
    output(:,1) = continuousInterpolant.X(position(:,2),position(:,1));
    output(:,2) = continuousInterpolant.Y(position(:,2),position(:,1));
else
    output(:,1) = etaInterpolant.X(position(:,2),position(:,1));
    output(:,2) = etaInterpolant.Y(position(:,2),position(:,1));
end

output = transpose(output);

% Orientation discontinuity
if ~isempty(previousEta.value) && ~all(isnan(previousEta.value))
    output = sign(previousEta.value*output)*output;
end

function status = ode_output(~,position,flag,eta,etaInterpolant,domain,...
    flowResolution,etaGrid)

if nargin < 3 || isempty(flag)

    % Use the last position
    position = position(:,end);
    
    continuousInterpolant = is_element_with_orient_discont(position,...
        domain,flowResolution,etaGrid);

    if ~isempty(continuousInterpolant)
        currentEta = [continuousInterpolant.X(position(2),position(1)) ...
            continuousInterpolant.Y(position(2),position(1))];
    else
        currentEta = [etaInterpolant.X(position(2),position(1)) ...
            etaInterpolant.Y(position(2),position(1))];
    end
    
    eta.value = sign(currentEta*transpose(eta.value))*currentEta;

else
    switch(flag)
        case 'init'
            
            continuousInterpolant = is_element_with_orient_discont(...
                position,domain,flowResolution,etaGrid);
            if ~isempty(continuousInterpolant)
                eta.value = [continuousInterpolant.X(position(2),...
                    position(1)) continuousInterpolant.Y(position(2),...
                    position(1))];
            else
                eta.value = [etaInterpolant.X(position(2),position(1)) ...
                    etaInterpolant.Y(position(2),position(1))];
            end
            
        case 'done'
            
    end
end

status = 0;

function [distance,isTerminal,direction] = ode_events(~,position,domain)

isTerminal = true;
direction = 1;

if any(isnan(position))
    distance = 0;
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
    
    continuousInterpolant.X = griddedInterpolant({positionY,positionX},...
        [vector3(1) vector4(1); vector2(1) vector1(1)]);
    continuousInterpolant.Y = griddedInterpolant({positionY,positionX},...
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
