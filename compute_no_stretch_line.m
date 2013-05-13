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
        warning([mfilename,':alpha_not_real'],'alpha not real')
    end
    
    beta = sqrt((1-l1)./(l2-l1));
    if ~isreal(beta)
        warning([mfilename,':beta_not_real'],'beta not real')
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

if ~isfield(noStretchLine,'initialPosition')
    noStretchLine.initialPosition = initialize_ic_grid(noStretchLine.resolution,flow.domain);
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

% Vector field must be normalized
finalTime = noStretchLine.maxLength;
timespan = [0 finalTime];

spmdFun = @(idx)integrate_line(timespan,...
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

spmdFun = @(idx)integrate_line(timespan,...
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

spmdFun = @(idx)integrate_line(timespan,...
    noStretchLine.initialPosition(idx,:),flow.domain,flow.resolution,...
    noStretchLine.chiPos,noStretchLine.odeSolverOptions);

spmd
    iNoStretchLineC = codistributed(iNoStretchLine);
    iNoStretchLineL = getLocalPart(iNoStretchLineC);
    positionPosL = arrayfun(spmdFun,iNoStretchLineL,'UniformOutput',false);
    positionPosC = codistributed.build(positionPosL,...
        getCodistributor(iNoStretchLineC));
end

spmdFun = @(idx)integrate_line(timespan,...
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
