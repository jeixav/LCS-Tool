function cgStrain_convergence

%% Input parameters
epsilon = .1;
amplitude = .1;
omega = pi/5;
domain = [0,2;0,1];
resolutionX = 2000:500:4000;
timespan = [0,5];

%% Velocity definition
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
odeSolverOptions = odeset('relTol',1e-3,'absTol',1e-6);
eigenvalueFromMainGrid = true;

for m = 1:numel(resolutionX)
    % Make x and y grid spacing as equal as possible
    lResolutionX = resolutionX(m);
    gridSpace = diff(domain(1,:))/(double(lResolutionX) - 1);
    resolutionY = round(diff(domain(2,:))/gridSpace) + 1;
    resolution = [lResolutionX,resolutionY];
    disp(['Resolution: ',num2str(resolution)])
    
    %% Cauchy-Green strain eigenvalues and eigenvectors
    cgEigenvalue = eig_cgStrain(lDerivative,domain,resolution,timespan,'odeSolverOptions',odeSolverOptions,'eigenvalueFromMainGrid',eigenvalueFromMainGrid);
    prodCgEigenvalue = prod(cgEigenvalue,2);
    fprintf('prod(cgEigenvalue): min = %g, max = %g\n',min(prodCgEigenvalue),max(prodCgEigenvalue))
%     hAxes = setup_figure(domain);
%     hProdCgEigenvalue = imagesc(domain(1,:),domain(2,:),reshape(prodCgEigenvalue,fliplr(resolution)));
%     set(hProdCgEigenvalue,'Parent',hAxes)
%     hColorbar = colorbar('peer',hAxes);
%     set(get(hColorbar,'xlabel'),'string','\lambda_1 \lambda_2')
%     drawnow
end
