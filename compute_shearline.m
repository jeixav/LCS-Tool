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
    progressBar.setElapsedTimeVisible(1)
    progressBar.setRemainedTimeVisible(1)
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

spmdFunPos = @(idx)integrate_line(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaPos,shearline.odeSolverOptions);
spmdFunNeg = @(idx)integrate_line(timespan,...
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
spmdFunPos = @(idx)integrate_line(timespan,...
    shearline.initialPosition(idx,:),flow.domain,flow.resolution,...
    shearline.etaPos,shearline.odeSolverOptions);
spmdFunNeg = @(idx)integrate_line(timespan,...
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
