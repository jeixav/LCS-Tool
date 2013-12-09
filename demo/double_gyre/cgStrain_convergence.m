function cgStrain_convergence

%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolutionX = 2000:250:2500;
timespan = [0,6];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
incompressible = false;

%% LCS parameters
cgStrainOdeSolverOptions = odeset('relTol',1e-3);
eigenvalueFromMainGrid = true;

for m = 1:numel(resolutionX)
    % Make x and y grid spacing as equal as possible
    lResolutionX = resolutionX(m);
    gridSpace = diff(domain(1,:))/(double(lResolutionX) - 1);
    resolutionY = round(diff(domain(2,:))/gridSpace) + 1;
    resolution = [lResolutionX,resolutionY];
    disp(['Resolution: ',num2str(resolution)])
    
    %% Cauchy-Green strain eigenvalues and eigenvectors
    cgEigenvalue = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions,'eigenvalueFromMainGrid',eigenvalueFromMainGrid);
    prodCgEigenvalue = prod(cgEigenvalue,2);
    fprintf('prod(cgEigenvalue): min = %g, max = %g\n',min(prodCgEigenvalue),max(prodCgEigenvalue))
end
