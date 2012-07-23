function strainline = compute_strainline(flow,strainline)

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

spmd
    iStrainlineC = codistributed(iStrainline);
    iStrainlineL = getLocalPart(iStrainlineC);
    positionL = arrayfun(spmdFun,iStrainlineL,'UniformOutput',false);
    positionC = codistributed.build(positionL,...
        getCodistributor(iStrainlineC));
end

strainline.position = gather(positionC);

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

% Concatenate forward and backward time integration results
tmp = cellfun(@(input)flipud(input(2:end,:)),gather(positionC),...
    'uniformOutput',false);
strainline.position = cellfun(@(a,b)[a;b],tmp,strainline.position,...
    'UniformOutput',false);
