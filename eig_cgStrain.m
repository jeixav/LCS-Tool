% Calculate eigenvalues and eigenvectors of Cauchy-Green strain

function [cgStrainD,cgStrainV,cgStrain,finalPosition,dFlowMap] = eig_cgStrain(flow,method,eigMethod,coupledIntegration,verbose)

narginchk(1,5)

if nargin < 2
    method.name = 'equationOfVariation';
end

if nargin < 3
    eigMethod = 'standard';
end

if nargin < 4
    coupledIntegration = 'false';
end

if nargin < 5
    verbose.progress = false;
    verbose.stats = false;
end

initialPosition = initialize_ic_grid(flow.resolution,flow.domain);

switch method.name
    case 'finiteDifference'
        
        if ~isfield(method,'auxiliaryGridRelativeDelta') || ...
                isempty(method.auxiliaryGridRelativeDelta)
            method.auxiliaryGridRelativeDelta = 1e-2;
            warning([mfilename,':defaultAuxiliaryGridRelativeDelta'],...
                ['auxiliaryGridRelativeDelta not set; using default value: ',...
                num2str(method.auxiliaryGridRelativeDelta)])
        end
        
        % Eigenvectors from auxiliary grid
        deltaX = (flow.domain(1,2) - flow.domain(1,1))...
            /double(flow.resolution(1))*method.auxiliaryGridRelativeDelta;
        auxiliaryGridAbsoluteDelta = deltaX;
        auxiliaryPosition = auxiliary_position(initialPosition,...
            auxiliaryGridAbsoluteDelta);
        
        % Transform auxiliaryPosition into a two column array
        auxiliaryPositionX = auxiliaryPosition(:,1:2:end-1);
        auxiliaryPositionY = auxiliaryPosition(:,2:2:end);
        auxiliaryPosition = [auxiliaryPositionX(:) auxiliaryPositionY(:)];
        
        finalPositionAuxGridSol = integrate_flow(flow,auxiliaryPosition);
        finalPositionAuxGrid = arrayfun(@(odeSolution)deval(odeSolution,...
            flow.timespan(2)),finalPositionAuxGridSol,'uniformOutput',...
            false);
        finalPositionAuxGrid = cell2mat(finalPositionAuxGrid);
        finalPositionAuxGrid = transpose(finalPositionAuxGrid);
        
        % Transform finalPosition into an eight column array
        finalPositionAuxGridX = finalPositionAuxGrid(:,1);
        finalPositionAuxGridY = finalPositionAuxGrid(:,2);
        nPoints = prod(double(flow.resolution));
        finalPositionAuxGridX = reshape(finalPositionAuxGridX,nPoints,4);
        finalPositionAuxGridY = reshape(finalPositionAuxGridY,nPoints,4);
        finalPositionAuxGrid = nan(nPoints,8);
        finalPositionAuxGrid(:,1:2:7) = finalPositionAuxGridX;
        finalPositionAuxGrid(:,2:2:8) = finalPositionAuxGridY;

        cgStrainAuxGrid = compute_cgStrain(finalPositionAuxGrid,flow,...
            method.auxiliaryGridRelativeDelta);
        
        [cgStrainV,cgStrainD] = arrayfun(@eig_array,...
            cgStrainAuxGrid(:,1),cgStrainAuxGrid(:,2),...
            cgStrainAuxGrid(:,3),'UniformOutput',false);
        
        cgStrainV = cell2mat(cgStrainV);
        
        if ~isfield(method,'eigenvalueFromMainGrid')
            method.eigenvalueFromMainGrid = true;
            warning([mfilename,':defaultEigenvalueFromMainGrid'],...
                ['eigenvalueFromMainGrid not set; using default value: ' ,...
                num2str(method.eigenvalueFromMainGrid)])
        end
        
        if method.eigenvalueFromMainGrid
            initialPosition = initialize_ic_grid(flow.resolution,flow.domain);
            
            finalPositionMainGridSol = integrate_flow(flow,initialPosition);
            finalPositionMainGrid = arrayfun(@(odeSolution)deval(odeSolution,...
            flow.timespan(2)),finalPositionMainGridSol,'uniformOutput',...
            false);
            finalPositionMainGrid = cell2mat(finalPositionMainGrid);
            finalPositionMainGrid = transpose(finalPositionMainGrid);
            
            cgStrainMainGrid = compute_cgStrain(finalPositionMainGrid,flow);
            
            [~,cgStrainD] = arrayfun(@eig_array,...
                cgStrainMainGrid(:,1),cgStrainMainGrid(:,2),...
                cgStrainMainGrid(:,3),'UniformOutput',false);
        end
        
        cgStrainD = cell2mat(cgStrainD);
        
        % Use the Cauchy-Green strain calculated with the auxiliary
        % grid for statistics.
        nRows = size(cgStrainAuxGrid,1);
        cgStrain = arrayfun(@(idx)[cgStrainAuxGrid(idx,1)...
            cgStrainAuxGrid(idx,2); cgStrainAuxGrid(idx,2)...
            cgStrainAuxGrid(idx,3)],1:nRows,'uniformOutput',false);
        cgStrain = cell2mat(cgStrain);
        cgStrain = reshape(cgStrain,[2 2 nRows]);

    case 'equationOfVariation'

        dFlowMap0 = eye(2);
        dFlowMap0 = reshape(dFlowMap0,4,1);
        nPosition = size(initialPosition,1);
        finalPosition = nan(nPosition,2);
        dFlowMap = nan(nPosition,4);

        if isfield(flow,'symDerivative') && ~isfield(flow,'dDerivative')
            symJacDy = symJacDerivative(flow.symDerivative);
            
            symVars = {'t','x','y'};
            jacDyScalar11 = matlabFunction(symJacDy{1,1},'vars',symVars);
            jacDyScalar12 = matlabFunction(symJacDy{1,2},'vars',symVars);
            jacDyScalar21 = matlabFunction(symJacDy{2,1},'vars',symVars);
            jacDyScalar22 = matlabFunction(symJacDy{2,2},'vars',symVars);
            
            flow.dDerivative = @(t,y)[jacDyScalar11(t,y(1),y(2)) ...
                jacDyScalar12(t,y(1),y(2)); jacDyScalar21(t,y(1),y(2)) ...
                jacDyScalar22(t,y(1),y(2))];
        end
        
        if ~isfield(flow,'odeSolver')
            odeSolver = @ode45;
            warning([mfilename,':defaultOdeSolver'],...
                ['odeSolver not set; using default: ',func2str(odeSolver)])
        else
            odeSolver = flow.odeSolver;
        end

        if isfield(flow,'odeSolverOptions')
            odeSolverOptions = flow.odeSolverOptions;
        else
            odeSolverOptions = [];
        end
        
        if ~isfield(flow,'derivative') && isfield(flow,'symDerivative')
            flow.derivative = sym2fun(flow.symDerivative);
        end
        
        if verbose.progress
            if ~exist('ParforProgressStarter2','file')
                addpath('ParforProgress2')
            end
            parforVerbose = true;
        else
            parforVerbose = false;
            progressBar = [];
        end
        
        if coupledIntegration
            % Add dFlowMap0
            initialPosition = [initialPosition,...
                repmat(transpose(dFlowMap0),size(initialPosition,1),1)];
            
            % Reshape array for pseudocoupled form
            initialPosition = transpose(initialPosition);
            initialPosition = initialPosition(:);
            
            % targetBlockSize controls total memory use; needs to be tuned
            % for different computers
            targetBlockSize = 50000;
            blockIndex = block_index(size(initialPosition,1),...
                targetBlockSize);
            
            nBlock = size(blockIndex,2);
            sol = cell(nBlock,1);

            if verbose.progress
                progressBar = ParforProgressStarter2(mfilename,nBlock);
            else
                progressBar = [];
            end
            
            parfor iBlock = 1:nBlock
                iBlockIndex = blockIndex(1,iBlock):blockIndex(2,iBlock); %#ok<PFBNS>
                [~,sol{iBlock}] = ode45(@(t,x)flow.derivative(t,x),flow.timespan,initialPosition(iBlockIndex),odeSolverOptions); %#ok<PFBNS>
                sol{iBlock} = sol{iBlock}(end,:);
                if ~isempty(progressBar)
                    progressBar.increment(iBlock) 
                end
            end
            sol = [sol{:}];
            
            sol = transpose(reshape(sol(end,:),6,size(sol,2)/6));
            dFlowMap = sol(:,3:6);
            % FIXME Check indices in flow definition file.
            dFlowMap(:,[2,3]) = fliplr(dFlowMap(:,[2,3]));
            finalPosition = sol(:,1:2);
        else
            if parforVerbose
                progressBar = ParforProgressStarter2(mfilename,...
                    nPosition);
            end
            parfor iPosition = 1:nPosition
                position0 = transpose(initialPosition(iPosition,:));
                y0 = [position0; dFlowMap0];
                sol = feval(odeSolver,@(t,y)eov_odefun(t,y,flow),...
                    flow.timespan,y0,odeSolverOptions);
                finalPosition(iPosition,:) = ...
                    transpose(deval(sol,flow.timespan(end),1:2));
                dFlowMap(iPosition,:) = ...
                    transpose(deval(sol,flow.timespan(end),3:6));
                if parforVerbose
                    progressBar.increment(iPosition) %#ok<PFBNS>
                end
            end
        end
                       
        if verbose.progress
            try
                delete(progressBar)
            catch me %#ok<NASGU>
            end
        end
                
        if isempty(eigMethod)
            eigMethod = 'standard';
            warning([mfilename,':defaultEigMethod'],...
                ['eigMethod not set; using default: ',eigMethod])
        end
        [cgStrainV,cgStrainD] = eov_compute_cgStrain(dFlowMap,eigMethod,...
            verbose);
        
        cgStrain = cgStrain_from_dFlowMap(dFlowMap);

end

if any(cgStrainD(:) <= 0)
    warning([mfilename,':nonpositiveEigenvalue'],...
        'Nonpositive eigenvalues')
end

if verbose.stats
    disp('cgStrain_stats:')
    cgStrain_stats(cgStrain,cgStrainV,cgStrainD,verbose.stats);
end

if ~isfield(flow,'imposeIncompressibility')
    flow.imposeIncompressibility = false;
    warning([mfilename,':imposeIncompressibility'],...
        ['imposeIncompressibility not set; using default value: ',...
        num2str(flow.imposeIncompressibility)])
end

if flow.imposeIncompressibility
    prodCgStrainD = prod(cgStrainD,2);
    if any(prodCgStrainD ~= 1)
        warning([mfilename,':eigenvalueProdNot1'],...
            'Eigenvalue products not 1')
        % Enforce incompressibility condition in eigenvalues
        cgStrainD(:,1) = 1./cgStrainD(:,2);
    end
end

function [v,d] = eig_array(x11,x12,x22)

[v,d] = eig([x11 x12; x12 x22]);

d = transpose(diag(d));
v = reshape(v,1,4);

function f = eov_odefun(t,y,flow)

f = nan(size(y));

f(1:2) = flow.derivative(t,y(1:2));

dFlowMap = reshape(y(3:6),2,2);
f2 = flow.dDerivative(t,y(1:2))*dFlowMap;
f(3:6) = reshape(f2,4,1);

% eov_compute_cgStrain Compute Cauchy-Green strain tensor, its eigenvalues
% and eigenvectors
%
% DESCRIPTION
% [cgStrainV,cgStrainD,cgStrain] = ...
%     eov_compute_cgStrain(dFlowMap,method,verbose)
%
% method is either 'standard' or 'custom'. 'standard' uses MATLAB's EIG
% function to calculate eigenvectors and eigenvalues. 'custom' calculates
% lambda_2 analytically, then lambda_1 = inv(lambda_2), then xi_2
% analytically and finally xi_1 = Omega*xi_2. This is only valid if the
% flow is incompressible.

function [cgStrainV,cgStrainD,cgStrain] = eov_compute_cgStrain(dFlowMap,...
    method,verbose)

nPosition = size(dFlowMap,1);
cgStrainV = nan(nPosition,4);
cgStrainD = nan(nPosition,2);
cgStrain = cell(nPosition,1);

if verbose.progress
    progressBar = ParforProgressStarter2(mfilename,nPosition);
    parforVerbose = true;
else
    parforVerbose = false;
    progressBar = [];
end

switch method
    case 'standard'
        parfor i = 1:nPosition
            dFlowMap2 = reshape(dFlowMap(i,:),2,2);
            cgStrain{i} = transpose(dFlowMap2)*dFlowMap2;
            [v,d] = eig(cgStrain{i});
            cgStrainV(i,:) = reshape(v,1,4);
            cgStrainD(i,:) = [d(1) d(4)];
            if parforVerbose
                progressBar.increment(i) %#ok<PFBNS>
            end
        end
    case 'custom'
        parfor i = 1:nPosition
            dFlowMap2 = reshape(dFlowMap(i,:),2,2);
            cgStrain{i} = transpose(dFlowMap2)*dFlowMap2;
            [v,d] = eig_custom(cgStrain{i});
            cgStrainV(i,:) = reshape(v,1,4);
            cgStrainD(i,:) = [d(1) d(4)];
            if parforVerbose
                progressBar.increment(i) %#ok<PFBNS>
            end
        end
    otherwise
        error('Invalid method')
end

% parfor i = 1:nPosition
%     dFlowMap2 = reshape(dFlowMap(i,:),2,2);
%     cgStrain{i} = transpose(dFlowMap2)*dFlowMap2;
%     switch method
%         case 'standard'
%             [v,d] = eig(cgStrain{i});
%         case 'custom'
%             [v,d] = eig_custom(cgStrain{i});
%         otherwise
%             error('Invalid method')
%     end
%     cgStrainV(i,:) = reshape(v,1,4);
%     cgStrainD(i,:) = [d(1) d(4)];
%     if parforVerbose
%         progressBar.increment(i) %#ok<PFBNS>
%     end
% end

if verbose.progress
    try
        delete(progressBar);
    catch me %#ok<NASGU>
    end
end

function cgStrain = cgStrain_from_dFlowMap(dFlowMap)

nRows = size(dFlowMap,1);
dFlowMap = reshape(transpose(dFlowMap),[2 2 nRows]);
cgStrain = arrayfun(@(idx)transpose(dFlowMap(:,:,idx))...
    *dFlowMap(:,:,idx),1:nRows,'UniformOutput',false);
cgStrain = cell2mat(cgStrain);
cgStrain = reshape(cgStrain,[2 2 nRows]);

function [v,d] = eig_custom(a)

d(2,2) = .5*trace(a) + sqrt((.5*trace(a)).^2 - det(a));
if imag(d(2,2))
    warning([mfilename,':complexEigenvalue'],['Complex eigenvalue: ',num2str(d(2,2))])
end
d(1,1) = inv(d(2,2));

if d(2,2) < d(1,1)
    warning([mfilename,':eigenvalueOrder','Eigenvalue ordering error'])
end

v(1,2) = -a(1,2)/sqrt(a(1,2)^2 + (a(1,1) - d(2,2))^2);
v(2,2) = (a(1,1) - d(2,2))/sqrt(a(1,2)^2 + (a(1,1) - d(2,2))^2);

v(1,1) = v(2,2);
v(2,1) = -v(1,2);

% Calculate block indices to perform for-loop integration
function blockIndex = block_index(nInitialPosition,targetBlockSize)

if mod(nInitialPosition,6)
    error('nInitialPosition')
end

blockSize = targetBlockSize - rem(targetBlockSize,6);

blockStartIndex = 1:blockSize:nInitialPosition;
blockEndIndex = [blockStartIndex(2:end)-1 nInitialPosition];

blockIndex = [blockStartIndex; blockEndIndex];
