function [cgStrainD,cgStrainV] = eig_cgStrain(flow,method,verbose)
% Calculate eigenvalues and eigenvectors of Cauchy-Green strain
%
% eig_cgStrain(flow)
% eig_cgStrain(flow,method)
% eig_cgStrain(flow,method,verbose)
%
% method is a structure. method.name can be 'fd' or 'eov'. If it is 'fd',
% method.eigenvalueFromMainGrid can be true or false.

narginchk(1,3)

if nargin < 3
    verbose.progress = false;
    verbose.stats = false;
end

if nargin < 2
    method.name = 'eov';
end

initialPosition = initialize_ic_grid(flow.resolution,flow.domain);

switch method.name
    case 'fd'
        
        % Eigenvectors from auxiliary grid
        deltaX = (flow.domain(1,2) - flow.domain(1,1))...
            /double(flow.resolution(1))*flow.auxiliaryGridRelativeDelta;
        delta = deltaX;
        auxiliaryPosition = auxiliary_position(initialPosition,delta);
        
        if verbose.progress
            fprintf('Auxilary grid\n')
            flow.odeSolverOptions.OutputFcn = @ode_progress_bar;
        end
        finalPositionAuxGrid = integrate_flow(flow,auxiliaryPosition);
        
        cgStrainAuxGrid = compute_cgStrain(finalPositionAuxGrid,flow);
        
        [cgStrainV,~] = arrayfun(@eig_array,...
            cgStrainAuxGrid(:,1),cgStrainAuxGrid(:,2),cgStrainAuxGrid(:,3),...
            'UniformOutput',false);
        
        cgStrainV = cell2mat(cgStrainV);
        
        if method.params.eigenvalueFromMainGrid
            % Eigenvalues from main grid
            initialPosition = initialize_ic_grid(flow.resolution,flow.domain);
            
            if verbose.progress
                fprintf('Main grid\n')
                flow.odeSolverOptions.OutputFcn = @ode_progress_bar;
            end
            finalPositionMainGrid = integrate_flow(flow,initialPosition);
            
            cgStrainMainGrid = compute_cgStrain(finalPositionMainGrid,flow);
            
            [~,cgStrainD] = arrayfun(@eig_array,...
                cgStrainMainGrid(:,1),cgStrainMainGrid(:,2),...
                cgStrainMainGrid(:,3),'UniformOutput',false);
        end
        
        cgStrainD = cell2mat(cgStrainD);
        
    case 'eov'

        dFlowMap0 = eye(2);
        dFlowMap0 = reshape(dFlowMap0,4,1);
        nPosition = size(initialPosition,1);
        dFlowMap = nan(nPosition,4);

        if isfield(flow,'odeSolverOptions')
            odeSolverOptions = flow.odeSolverOptions;
        else
            odeSolverOptions = [];
        end
        odeSolver = flow.odeSolver;

        if verbose.progress
            progressBar = ConsoleProgressBar;
            progressBar.setText([mfilename,' ODE'])
            progressBar.setTextPosition('left')
            progressBar.setElapsedTimeVisible(1)
            progressBar.setRemainedTimeVisible(1)
            progressBar.setLength(20)
            progressBar.setMaximum(nPosition)
            progressBar.start
        end

        parfor iPosition = 1:nPosition
            position0 = transpose(initialPosition(iPosition,:));
            y0 = [position0; dFlowMap0];
            sol = feval(odeSolver,@(t,y)eov_odefun(t,y,flow),...
                flow.timespan,y0,odeSolverOptions);
            dFlowMap(iPosition,:) = ...
                transpose(deval(sol,flow.timespan(end),3:6));
            if verbose.progress && matlabpool('size') == 0
                progressBar.setValue(progressBar.value + 1)
            end
        end

        if verbose.progress
            progressBar.setValue(progressBar.maximum)
            progressBar.stop
            fprintf('\n')
        end
        
        [cgStrainV,cgStrainD,cgStrain] = eov_compute_cgStrain(dFlowMap,...
            verbose);
end

if any(cgStrainD(:) <= 0)
    warning('eig_cgStrain:nonpositiveEigenvalue','Nonpositive eigenvalues')
end

if ~flow.isCompressible
    prodCgStrainD = prod(cgStrainD,2);
    if any(prodCgStrainD ~= 1)
        warning('eig_cgStrain:eigenvalueProdNot1',...
            'Eigenvalue products not 1')
        if verbose.stats
            fprintf('min = %g\n',min(prodCgStrainD))
            fprintf('max = %g\n',max(prodCgStrainD))
            fprintf('mean = %g\n',mean(prodCgStrainD))
            fprintf('median = %g\n',median(prodCgStrainD))
            fprintf('\n')
            fprintf('max lambda_1 = %g\n',max(cgStrainD(:,1)))
            detCgStrain = cellfun(@det,cgStrain);
            fprintf('mean(abs(detCgStrain-1)) = %g\n',...
                mean(abs(detCgStrain-1)))
        end
        
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

function [cgStrainV,cgStrainD,cgStrain] = eov_compute_cgStrain(dFlowMap,...
    verbose)

nPosition = size(dFlowMap,1);
cgStrainV = nan(nPosition,4);
cgStrainD = nan(nPosition,2);
cgStrain = cell(nPosition,1);

if verbose.progress
    progressBar = ConsoleProgressBar;
    progressBar.setText([mfilename,' EIG'])
    progressBar.setTextPosition('left')
    progressBar.setElapsedTimeVisible(1)
    progressBar.setRemainedTimeVisible(1)
    progressBar.setLength(20)
    progressBar.setMaximum(nPosition)
    progressBar.start
end

parfor i = 1:nPosition
    dFlowMap2 = reshape(dFlowMap(i,:),2,2);
    cgStrain{i} = transpose(dFlowMap2)*dFlowMap2;
    [v,d] = eig(cgStrain{i});
    cgStrainV(i,:) = reshape(v,1,4);
    cgStrainD(i,:) = [d(1) d(4)];
    if verbose.progress && matlabpool('size') == 0
        progressBar.setValue(progressBar.value + 1)
    end
end

if verbose.progress
    progressBar.setValue(progressBar.maximum)
    progressBar.stop
    fprintf('\n')
end
