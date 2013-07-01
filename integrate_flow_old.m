%integrate_flow Integrate flow
%
% SYNTAX
% flowSolution = integrate_flow(flow,initialPosition,useEoV)
% flowSolution = integrate_flow(flow,initialPosition,useEoV,verbose)
%
% INPUT ARGUMENTS
% initialPosition: n-by-2 array
% useEoV: true or false
% verbose: true or false

function flowSolution = integrate_flow(flow,initialPosition,useEoV,verbose)

if nargin < 4
    verbose = true;
end

nPosition = size(initialPosition,1);

if verbose
    if ~exist('ParforProgressStarter2','file')
        addpath('ParforProgress2')
    end
    progressBar = ParforProgressStarter2(mfilename,nPosition);
else
    progressBar = [];
end

if ~isfield(flow,'odeSolverOptions')
    odeSolverOptions = [];
else
    odeSolverOptions = flow.odeSolverOptions;
end
odefun = flow.derivative;
timespan = flow.timespan;

if flow.coupledIntegration
    % FIXME Copy-paste from eig_cgStrain
    initialPosition = transpose(initialPosition);
    initialPosition = initialPosition(:);
    
    targetBlockSize = 1e9;
    blockIndex = block_index(size(initialPosition,1),targetBlockSize);

    nBlock = size(blockIndex,2);
    flowSolution = cell(nBlock,1);
            
    ticID = tic;
    disp([mfilename,' progress:'])
    reverseStr = '';
    for iBlock = 1:nBlock
        iBlockIndex = blockIndex(1,iBlock):blockIndex(2,iBlock);
        flowSolution{iBlock} = ode45(@(t,y)odefun(t,y,useEoV),timespan,initialPosition(iBlockIndex),odeSolverOptions);
        elapsed = toc(ticID);
        total = toc(ticID)/(iBlock/nBlock);
        msg = sprintf('Elapsed: %s Remaing: %s Total: %s',seconds2human(elapsed,'full'),seconds2human(total-elapsed,'short'),seconds2human(total,'short'));
        fprintf([reverseStr,msg])
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    flowSolution = [flowSolution{:}];
else
    flowCgStrainMethodName = flow.cgStrainMethod.name;

    % FIXME If ODE solution structure is changed by MathWorks, this will
    % break
    flowSolution(nPosition) = struct('solver',[],'extdata',[],'x',[],'y',[],'stats',[],'idata',[]);

    switch flowCgStrainMethodName
        case 'finiteDifference'
            parfor iPosition = 1:nPosition
                flowSolution(iPosition) = ode45(@(t,y)odefun(t,y,useEoV),timespan,initialPosition(iPosition,:),odeSolverOptions); %#ok<PFBNS>
                if verbose
                    progressBar.increment(iPosition) %#ok<PFBNS>
                end
            end
        case 'equationOfVariation'
            parfor iPosition = 1:nPosition
                if useEoV
                    dFlowMap0 = eye(2);
                    dFlowMap0 = reshape(dFlowMap0,4,1);
                    iInitialPosition = [initialPosition(iPosition,:),dFlowMap0];
                else
                    iInitialPosition = initialPosition(iPosition,:);
                end
                flowSolution(iPosition) = ode45(@(t,y)odefun(t,y,useEoV),timespan,iInitialPosition,odeSolverOptions); %#ok<PFBNS>
                if verbose
                    progressBar.increment(iPosition) %#ok<PFBNS>
                end
            end
    end
end

if verbose
    try
        delete(progressBar)
    catch me %#ok<NASGU>
    end
end

% Calculate block indices to perform hybrid vector/for-loop integration
% FIXME Copy-paste from eig_cgStrain
function blockIndex = block_index(nInitialPosition,targetBlockSize)

if mod(nInitialPosition,2)
    error('nInitialPosition')
end

blockSize = targetBlockSize - rem(targetBlockSize,6);

blockStartIndex = 1:blockSize:nInitialPosition;
blockEndIndex = [blockStartIndex(2:end)-1 nInitialPosition];

blockIndex = [blockStartIndex; blockEndIndex];
