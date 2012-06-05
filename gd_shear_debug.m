function gd_shear_debug

addpath('flow_templates')

doubleGyre = double_gyre;

geodesicDeviationPosI = computation(doubleGyre);

doubleGyre.flow.isCompressible = true;
geodesicDeviationPosC = computation(doubleGyre);

positionX = linspace(doubleGyre.flow.domain(1,1),...
    doubleGyre.flow.domain(1,2),doubleGyre.flow.resolution(1));
positionY = linspace(doubleGyre.flow.domain(2,1),...
    doubleGyre.flow.domain(2,2),doubleGyre.flow.resolution(2));

figure
contourf(positionX,positionY,geodesicDeviationPosI,linspace(1e-8,1,20),...
    'lineStyle','none')
colorbar
set(gca,'DataAspectRatio',[1 1 1])
title('Incompressible')

fprintf('Incompressible:\n')
print_stats(geodesicDeviationPosI)

figure
contourf(positionX,positionY,geodesicDeviationPosC,linspace(1e-8,1,20),...
    'lineStyle','none')
colorbar
set(gca,'DataAspectRatio',[1 1 1])
title('Compressible')

fprintf('Compressible:\n')
print_stats(geodesicDeviationPosC)

figure
contourf(positionX,positionY,...
    abs(geodesicDeviationPosI-geodesicDeviationPosC),...
    linspace(5e-7,1,20),'lineStyle','none')
colorbar
set(gca,'DataAspectRatio',[1 1 1])
title('Absolute difference')

function print_stats(geodesicDeviation)

fprintf('min = %g\n',min(geodesicDeviation(:)))
fprintf('max = %g\n',max(geodesicDeviation(:)))
fprintf('mean = %g\n',mean(geodesicDeviation(:)))
fprintf('median = %g\n',median(geodesicDeviation(:)))

function [geodesicDeviationPos,geodesicDeviationNeg] = computation(input)

flow = set_flow_default(input.flow);

% Computations
if ~isfield(flow,'finalPosition')
    initialPosition = initial_position(flow.domain,flow.resolution);
    deltaX = (flow.domain(1,2) - flow.domain(1,1))...
        /double(flow.resolution(1))*flow.auxiliaryGridRelativeDelta;
    delta = deltaX;
    auxiliaryPosition = auxiliary_position(initialPosition,delta);
    
    flow.odeSolverOptions.OutputFcn = @ode_progress_bar;
    flow.finalPosition = integrate_flow(flow,auxiliaryPosition);
end

if ~isfield(flow,'cgStrain')
    flow.cgStrain = compute_cgStrain(flow.finalPosition,...
        delta);
end

if ~all(isfield(flow,{'cgEigenvector','cgEigenvector'}))
    [flow.cgEigenvector,flow.cgEigenvalue] = arrayfun(@eig_array,...
        flow.cgStrain(:,1),flow.cgStrain(:,2),flow.cgStrain(:,3),...
        'UniformOutput',false);
    flow.cgEigenvalue = cell2mat(flow.cgEigenvalue);
    flow.cgEigenvector = cell2mat(flow.cgEigenvector);
end

if ~all(isfield(shearline,{'etaPos','etaNeg','positionPos','positionNeg'}))
    shearline.odeSolverOptions = odeset(shearline.odeSolverOptions,...
        'outputFcn',@ode_progress_bar);
    shearline = compute_shearline(flow,shearline);
end

if ~all(isfield(shearline,{'geodesicDeviationPos','geodesicDeviationNeg'}))
    [~,geodesicDeviationPos,geodesicDeviationNeg] = ...
        geodesic_deviation_shearline(flow,shearline,false);
end
