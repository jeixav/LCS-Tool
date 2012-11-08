% Compare convergence of Cauchy-Green strain between equation of variation
% method and finite difference with auxiliary grid method.

function data = compare_fd_eov

if ~matlabpool('size')
    matlabpool('open')
end

pctRunOnAll javaaddpath('ParforProgress2')

if ~exist('double_gyre','file')
    addpath('flow_templates')
end

doubleGyre = double_gyre;
doubleGyre.flow = set_flow_resolution([2 1]*10,doubleGyre.flow);

errTol = logspace(-2,-12,8);
cgStrainDetFD = nan(size(errTol));
cgStrainDetEoV = nan(size(errTol));

verbose.progress = true;
verbose.stats = false;

doubleGyre.flow.derivative = sym2fun(doubleGyre.flow.symDerivative);
doubleGyre.flow.odeSolver = @ode45;
    
for idx = 1:length(errTol)

    doubleGyre.flow.odeSolverOptions = odeset('relTol',errTol(idx),...
        'absTol',errTol(idx));

    % Finite difference method
    cgStrainMethod.name = 'finiteDifference';
    cgStrainMethod.auxiliaryGridRelativeDelta = 1e-8;
    cgStrainMethod.eigenvalueFromMainGrid = true;
    
    [cgEigenvalue,~,cgStrain] = eig_cgStrain(doubleGyre.flow,...
        cgStrainMethod,verbose);

    o = cgStrain_stats(cgStrain,cgEigenvalue,false);
    cgStrainDetFD(idx) = o(2);
    
    % Equation of variation method
    cgStrainMethod.name = 'equationOfVariation';
    [cgEigenvalue,~,cgStrain] = eig_cgStrain(doubleGyre.flow,...
        cgStrainMethod,verbose);

    o = cgStrain_stats(cgStrain,cgEigenvalue,false);
    cgStrainDetEoV(idx) = o(2);

end

data = [errTol;cgStrainDetFD;cgStrainDetEoV];

f1 = figure;
a1 = axes;
set(a1,'parent',f1)
set(a1,'nextplot','add')
set(a1,'xscale','log')
set(a1,'yscale','log')
set(a1,'box','on')
set(a1,'xgrid','on')
set(a1,'ygrid','on')

xlabel('Tolerance')
ylabel('Error')

p1 = plot(errTol,cgStrainDetFD);
set(p1,'marker','o')
set(p1,'displayName','Finite difference')

p2 = plot(errTol,cgStrainDetEoV);
set(p2,'marker','square')
set(p2,'displayName','Equation of variation')

set([p1 p2],'color','k')
set([p1 p2],'markerSize',8)
set([p1 p2],'markerFaceColor','k')

hLegend = legend([p1 p2]);
set(hLegend,'location','NorthWest')

% doubleGyre = strain_lcs_script(doubleGyre);
% 
% set(findobj(gca,'tag','strainline'),'visible','on')
