function [cgStrainD,cgStrainV] = eig_cgStrain(flow,...
    eigenvalueFromMainGrid,verbose)
% Calculate eigenvalues and eigenvectors of Cauchy-Green strain

if nargin < 3
    verbose.progressBar = false;
    verbose.stats = false;
end

if nargin < 2
    eigenvalueFromMainGrid = true;
end

% Eigenvectors from auxiliary grid
initialPosition = initial_position(flow.domain,flow.resolution);
deltaX = (flow.domain(1,2) - flow.domain(1,1))...
    /double(flow.resolution(1))*flow.auxiliaryGridRelativeDelta;
delta = deltaX;
auxiliaryPosition = auxiliary_position(initialPosition,delta);

if verbose.progressBar
    fprintf('Auxilary grid\n')
    flow.odeSolverOptions.OutputFcn = @ode_progress_bar;
end
finalPositionAuxGrid = integrate_flow(flow,auxiliaryPosition);

cgStrainAuxGrid = compute_cgStrain(finalPositionAuxGrid,flow);

[cgStrainV,~] = arrayfun(@eig_array,...
    cgStrainAuxGrid(:,1),cgStrainAuxGrid(:,2),cgStrainAuxGrid(:,3),...
    'UniformOutput',false);

cgStrainV = cell2mat(cgStrainV);

if eigenvalueFromMainGrid
    % Eigenvalues from main grid
    initialPosition = initial_position(flow.domain,flow.resolution);
    
    if verbose.progressBar
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
        end
    end
end

end

function [v,d] = eig_array(x11,x12,x22)

[v,d] = eig([x11 x12; x12 x22]);

d = transpose(diag(d));
v = reshape(v,1,4);

end
