function shearline = compute_shearline(flow,shearline,verbose)

narginchk(2,3)

if nargin < 3
    verbose.progress = false;
end

if ~all(isfield(shearline,{'etaPos','etaNeg'}))
    [shearline.etaPos,shearline.etaNeg] = lagrangian_shear(flow.cgEigenvector,flow.cgEigenvalue);
    if ~isreal(shearline.etaPos)
        warning([mfilename,':etaPos_not_real'],'etaPos not real')
    end
    if ~isreal(shearline.etaNeg)
        warning([mfilename,':etaNeg_not_real'],'etaNeg not real')
    end
end

% Vector field must be normalized
finalTime = shearline.maxLength;
timespan = [0 finalTime];

if ~isfield(shearline,'initialPosition')
    shearline.initialPosition = initialize_ic_grid(shearline.resolution,flow.domain);
end

nShearlines = size(shearline.initialPosition,1);

if verbose.progress
    if ~exist('ParforProgressStarter2','file')
        addpath('ParforProgress2')
    end
    progressBar = ParforProgressStarter2(mfilename,2*nShearlines);
else
    progressBar = false;
end

if ~isfield(shearline,'odeSolverOptions')
    odeSolverOptions = odeset;
else
    odeSolverOptions = odeset(shearline.odeSolverOptions);
end

parforFunPos = @(idx)integrate_line(timespan,shearline.initialPosition(idx,:),flow.domain,flow.resolution,flow.periodicBc,shearline.etaPos,odeSolverOptions);
parforFunNeg = @(idx)integrate_line(timespan,shearline.initialPosition(idx,:),flow.domain,flow.resolution,flow.periodicBc,shearline.etaNeg,odeSolverOptions);

parfor i = 1:nShearlines
    positionPos{i} = feval(parforFunPos,i);
    positionNeg{i} = feval(parforFunNeg,i);
    if verbose.progress %#ok<PFBNS>
        progressBar.increment(i) %#ok<PFBNS>
    end
end

shearline.positionPos = positionPos;
shearline.positionNeg = positionNeg;

% Backward time integration
timespan = -timespan;
parforFunPos = @(idx)integrate_line(timespan,shearline.initialPosition(idx,:),flow.domain,flow.resolution,flow.periodicBc,shearline.etaPos,odeSolverOptions);
parforFunNeg = @(idx)integrate_line(timespan,shearline.initialPosition(idx,:),flow.domain,flow.resolution,flow.periodicBc,shearline.etaNeg,odeSolverOptions);

parfor i = 1:nShearlines
    positionPos{i} = feval(parforFunPos,i);
    positionNeg{i} = feval(parforFunNeg,i);
    if verbose.progress %#ok<PFBNS>
        progressBar.increment(i+nShearlines) %#ok<PFBNS>
    end
end

if verbose.progress
    try
        delete(progressBar);
    catch me %#ok<NASGU>
    end
end

% Concatenate forward and backward time integration results
tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionPos),...
    'uniformOutput',false);
shearline.positionPos = cellfun(@(a,b)[a;b],tmp,shearline.positionPos,...
    'UniformOutput',false);
tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionNeg),...
    'uniformOutput',false);
shearline.positionNeg = cellfun(@(a,b)[a;b],tmp,shearline.positionNeg,...
    'UniformOutput',false);
