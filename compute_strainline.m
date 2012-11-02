function strainline = compute_strainline(flow,strainline,verbose)

narginchk(2,3)

if nargin < 3
    verbose.progress = false;
end

% Vector field must be normalized
finalTime = strainline.maxLength;
timespan = [0 finalTime];

if isfield(strainline,'resolution')
    strainline.initialPosition = initialize_ic_grid(...
        strainline.resolution,flow.domain);
end

nStrainlines = size(strainline.initialPosition,1);

if ~isfield(strainline,'odeSolverOptions')
    strainline.odeSolverOptions = [];
end


if verbose.progress
    progressBar = ParforProgressStarter2(mfilename,2*nStrainlines);
end
        
parforFun = @(idx)integrate_line(timespan,...
    strainline.initialPosition(idx,:),flow.domain,flow.resolution,...
    flow.cgEigenvector(:,1:2),strainline.odeSolverOptions);

parfor i = 1:nStrainlines
    positionForward{i} = feval(parforFun,i);
    if verbose.progress %#ok<PFBNS>
        progressBar.increment(i) %#ok<PFBNS>
    end
end
        
strainline.position = positionForward;
        
timespan = -timespan;
parforFun = @(idx)integrate_line(timespan,...
    strainline.initialPosition(idx,:),flow.domain,flow.resolution,...
    flow.cgEigenvector(:,1:2),strainline.odeSolverOptions);

parfor i = 1:nStrainlines
    positionBackward{i} = feval(parforFun,i);
    if verbose.progress %#ok<PFBNS>
        progressBar.increment(i+nStrainlines) %#ok<PFBNS>
    end
end
        
if verbose.progress
    try
        delete(progressBar)
    catch me %#ok<NASGU>
    end
end

% Concatenate forward and backward time integration results
tmp = cellfun(@(input)flipud(input(2:end,:)),positionBackward,...
    'uniformOutput',false);
strainline.position = cellfun(@(a,b)[a;b],tmp,strainline.position,...
    'UniformOutput',false);
     
