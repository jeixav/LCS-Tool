function shearline = compute_shearline(flow,shearline,verbose)

narginchk(2,3)

if nargin < 3
    verbose.progress = false;
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

% Vector field must be normalized
finalTime = shearline.maxLength;
timespan = [0 finalTime];

if isfield(shearline,'resolution')
    shearline.initialPosition = initialize_ic_grid(shearline.resolution,...
        flow.domain);
end

nShearlines = size(shearline.initialPosition,1);

if verbose.progress
    progressBar = ParforProgressStarter2(mfilename,2*nShearlines);
else
    progressBar = false;
end

if ~isfield(shearline,'odeSolverOptions')
    shearline.odeSolverOptions = [];
end

parforFunPos = @(idx)integrate_line(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaPos,shearline.odeSolverOptions);
parforFunNeg = @(idx)integrate_line(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaNeg,shearline.odeSolverOptions);

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
parforFunPos = @(idx)integrate_line(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaPos,shearline.odeSolverOptions);
parforFunNeg = @(idx)integrate_line(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaNeg,shearline.odeSolverOptions);

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
