function strainline = compute_strainline(flow,strainline,verbose)

if nargin < 3
    verbose.progress = false;
end

timespan = [0 strainline.finalTime];

if isfield(strainline,'resolution')
    strainline.initialPosition = initialize_ic_grid(...
        strainline.resolution,flow.domain);
end

nStrainlines = size(strainline.initialPosition,1);
iStrainline = 1:nStrainlines;

if ~isfield(strainline,'odeSolverOptions')
    strainline.odeSolverOptions = [];
end

spmdFun = @(idx)integrate_line(timespan,...
    strainline.initialPosition(idx,:),flow.domain,flow.resolution,...
    flow.cgEigenvector(:,1:2),strainline.odeSolverOptions);

if verbose.progress
    progressBar = ConsoleProgressBar;
    progressBar.setText(mfilename)
    progressBar.setTextPosition('left')
    progressBar.setElapsedTimeVisible(1)
    progressBar.setRemainedTimeVisible(1)
    progressBar.setLength(20)
    progressBar.setMaximum(2*nStrainlines)
    progressBar.start
end

spmd
    iStrainlineC = codistributed(iStrainline);
    iStrainlineL = getLocalPart(iStrainlineC);
    positionL = arrayfun(spmdFun,iStrainlineL,'UniformOutput',false);
    positionC = codistributed.build(positionL,...
        getCodistributor(iStrainlineC));
end

strainline.position = gather(positionC);
progressBar.setValue(nStrainlines)

% Backward time integration
timespan = -timespan;

spmdFun = @(idx)integrate_line(timespan,...
    strainline.initialPosition(idx,:),flow.domain,flow.resolution,...
    flow.cgEigenvector(:,1:2),strainline.odeSolverOptions);

spmd
    iStrainlineC = codistributed(iStrainline);
    iStrainlineL = getLocalPart(iStrainlineC);
    positionL = arrayfun(spmdFun,iStrainlineL,'UniformOutput',false);
    positionC = codistributed.build(positionL,...
        getCodistributor(iStrainlineC));
end

if verbose.progress
    progressBar.setValue(progressBar.maximum)
    progressBar.stop
    fprintf('\n')
end
        
% Concatenate forward and backward time integration results
tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionC),...
    'uniformOutput',false);
strainline.position = cellfun(@(a,b)[a;b],tmp,strainline.position,...
    'UniformOutput',false);
