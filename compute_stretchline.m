function position = compute_stretchline(flow,stretchline,varargin)

narginchk(2,3)

p = inputParser;
addOptional(p,'verbose',false,@islogical)

parse(p,varargin{:})
verbose = p.Results.verbose;

% Vector field must be normalized
finalTime = stretchline.maxLength;
timespan = [0 finalTime];

if isfield(stretchline,'initialPosition')
    initialPosition = stretchline.initialPosition;
else
    % Place initial positions within domain
    delta = diff(flow.domain,1,2)./transpose(double(stretchline.resolution) + 1);
    gridPosition{1} = linspace(flow.domain(1,1) + .5*delta(1),flow.domain(1,2) - .5*delta(1),stretchline.resolution(1));
    gridPosition{2} = linspace(flow.domain(2,1) + .5*delta(2),flow.domain(2,2) - .5*delta(2),stretchline.resolution(2));
    [gridPosition{1},gridPosition{2}] = meshgrid(gridPosition{1},gridPosition{2});
    initialPosition = [gridPosition{1}(:),gridPosition{2}(:)];
end

nStretchlines = size(initialPosition,1);

if ~isfield(stretchline,'odeSolverOptions')
    odeSolverOptions = odeset;
else
    odeSolverOptions = odeset(stretchline.odeSolverOptions);
end

if verbose
    if ~exist('ParforProgressStarter2','file')
        addpath('ParforProgress2')
    end
    progressBar = ParforProgressStarter2(mfilename,2*nStretchlines);
end

flowPeriodicBc = [false,false];
parforFun = @(idx)integrate_line(timespan,initialPosition(idx,:),flow.domain,flow.resolution,flowPeriodicBc,flow.cgEigenvector(:,3:4),odeSolverOptions);

if ~verbose
    progressBar = [];
end

parfor i = 1:nStretchlines
    positionForward{i} = feval(parforFun,i);
    if verbose
        progressBar.increment(i) %#ok<PFBNS>
    end
end

position = positionForward;

timespan = -timespan;
parforFun = @(idx)integrate_line(timespan,initialPosition(idx,:),flow.domain,flow.resolution,flowPeriodicBc,flow.cgEigenvector(:,3:4),odeSolverOptions);

parfor i = 1:nStretchlines
    positionBackward{i} = feval(parforFun,i);
    if verbose
        progressBar.increment(i+nStretchlines) %#ok<PFBNS>
    end
end

if verbose
    try
        delete(progressBar)
    catch me %#ok<NASGU>
    end
end

% Concatenate forward and backward time integration results
tmp = cellfun(@(input)flipud(input(2:end,:)),positionBackward,'uniformOutput',false);
position = cellfun(@(a,b)[a;b],tmp,position,'uniformOutput',false);
