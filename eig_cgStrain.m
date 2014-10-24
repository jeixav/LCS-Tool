% eig_cgStrain Calculate eigenvalues and eigenvectors of Cauchy-Green strain
%
% SYNTAX
% cgStrainD = eig_cgStrain(derivative,domain,timespan,resolution)
% [cgStrainV,cgStrainD] = eig_cgStrain(derivative,domain,timespan,resolution)
% [cgStrainV,cgStrainD] = eig_cgStrain(...,'auxGridRelDelta',auxGridRelDelta)
% [cgStrainV,cgStrainD] = eig_cgStrain(...,'eigenvalueFromMainGrid',eigenvalueFromMainGrid)
% [cgStrainV,cgStrainD] = eig_cgStrain(...,'incompressible',icompressible)
% [cgStrainV,cgStrainD] = eig_cgStrain(...,'odeSolverOptions',options)
%
% INPUT ARGUMENTS
% derivative: function handle that evaluates the flow velocity equations.
% The function must have the form:
% [velX1;velY1;velX2;velY2;...;velXN;velYN] = f(t,[x1;y1;x2;y2;...;xN;yN].
% domain: defines the flow domain. It is of the form:
% [xMin,xMax;yMin,yMax].
% timespan: defines the flow timespan. It is of the form: [tStart,tEnd].
% resolution: defines the Cauchy-Green strain resolution. It is of the
% form: [resolutionX,resolutionY].
% auxGridRelDelta: number between 0 and 0.5 specifying the auxiliary grid
% spacing relative to the main grid spacing. Default is 1e-2.
% eigenvalueFromMainGrid: logical to control whether eigenvalues of the
% Cauchy-Green strain are calculated from main grid or auxiliary grid
% points. Default is true.
% incompressibile: logical that specifies if incompressibility should be
% imposed. If true, the smaller eigenvalue is set equal to the inverse of
% the larger eigenvalue, except at grid points where the larger eigenvalue
% is less than 1. At those grid points, eigenvalues and eigenvectors are
% set to NaN. Default is false.
% odeSolverOptions: odeset structure to specify flow map integration
% options.
%
% OUTPUT ARGUMENTS
% cgStrainD: prod(resolution)-by-2 array of Cauchy-Green strain
% eigenvalues.
% cgStrainV: prod(resolution)-by-4 array of Cauchy-Green strain
% eigenvectors.

function varargout = eig_cgStrain(derivative,domain,resolution,timespan,varargin)

nargoutchk(1,2)

p = inputParser;
addRequired(p,'derivative',@(derivative)validateattributes(derivative,{'function_handle'},{'scalar'}))
addRequired(p,'domain',@(domain)validateattributes(domain,{'double'},{'size',[2,2],'real','finite'}))
addRequired(p,'resolution',@(resolution)validateattributes(resolution,{'double'},{'size',[1,2],'real','finite'}))
addRequired(p,'timespan',@(timespan)validateattributes(timespan,{'double'},{'size',[1,2],'real','finite'}))
addParameter(p,'auxGridRelDelta',1e-2,@(input)validateattributes(input,{'double'},{'scalar','>',0,'<',.5}))
addParameter(p,'eigenvalueFromMainGrid',true,@(input)validateattributes(input,{'logical'},{'scalar'}))
addParameter(p,'incompressible',false,@(input)validateattributes(input,{'logical'},{'scalar'}))
addParameter(p,'odeSolverOptions',odeset)
% Undocumented parameters
addParameter(p,'customEigMethod',false,@(input)validateattributes(input,{'logical'},{'scalar'}))
addParameter(p,'coupledIntegration',true,@(input)validateattributes(input,{'logical'},{'scalar'}))
addParameter(p,'method','finiteDifference',@isstr)

parse(p,derivative,domain,resolution,timespan,varargin{:})

auxGridRelDelta = p.Results.auxGridRelDelta;
eigenvalueFromMainGrid = p.Results.eigenvalueFromMainGrid;
incompressible = p.Results.incompressible;
odeSolverOptions = p.Results.odeSolverOptions;
customEigMethod = p.Results.customEigMethod;
coupledIntegration = p.Results.coupledIntegration;
method = p.Results.method;

initialPosition = initialize_ic_grid(resolution,domain);

switch method
    case 'finiteDifference'
        
        if nargout == 2 || eigenvalueFromMainGrid == false
            %% Eigenvectors from auxiliary grid
            initialPositionM = reshape(initialPosition,[fliplr(resolution),2]);
            deltaX = mean(diff(initialPositionM(1,:,1)));
            deltaY = mean(diff(initialPositionM(:,1,2)));
            if deltaX ~= deltaY
                warning([mfilename,':unequalDelta'],'Unequal grid spacing: (deltaX - deltaY)/min([deltaX,deltaY]) = %.3g. Using deltaX to set auxiliary grid spacing.',(deltaX - deltaY)/min([deltaX,deltaY]))
            end
            gridSpace = deltaX;
            auxiliaryGridAbsoluteDelta = gridSpace*auxGridRelDelta;
            auxiliaryPosition = auxiliary_position(initialPosition,auxiliaryGridAbsoluteDelta);
            
            % Transform auxiliaryPosition into a two column array
            auxiliaryPositionX = auxiliaryPosition(:,1:2:end-1);
            auxiliaryPositionY = auxiliaryPosition(:,2:2:end);
            auxiliaryPositionC = [auxiliaryPositionX(:),auxiliaryPositionY(:)];
            
            if coupledIntegration
                finalPositionAuxGrid = ode45_vector(@(t,y)derivative(t,y,false),timespan,auxiliaryPositionC,false,odeSolverOptions);
            else
                error('Uncoupled integration not programmed.')
            end
            
            % Transform finalPosition into an eight column array
            finalPositionAuxGridX = finalPositionAuxGrid(:,1);
            finalPositionAuxGridY = finalPositionAuxGrid(:,2);
            nPoints = prod(double(resolution));
            finalPositionAuxGridX = reshape(finalPositionAuxGridX,nPoints,4);
            finalPositionAuxGridY = reshape(finalPositionAuxGridY,nPoints,4);
            finalPositionAuxGrid = nan(nPoints,8);
            finalPositionAuxGrid(:,1:2:7) = finalPositionAuxGridX;
            finalPositionAuxGrid(:,2:2:8) = finalPositionAuxGridY;
            
            cgStrainAuxGrid = compute_cgStrain(finalPositionAuxGrid,auxiliaryPosition,resolution);
            
            [cgStrainV,cgStrainD] = arrayfun(@(x11,x12,x22)eig_array(x11,x12,x22,customEigMethod),cgStrainAuxGrid(:,1),cgStrainAuxGrid(:,2),cgStrainAuxGrid(:,3),'UniformOutput',false);
            
            cgStrainV = cell2mat(cgStrainV);
        end
        
        %% Eigenvalues from main grid
        if eigenvalueFromMainGrid
            if coupledIntegration
                finalPositionMainGrid = ode45_vector(@(t,y)derivative(t,y,false),timespan,initialPosition,false,odeSolverOptions);
            else
                finalPositionMainGrid = nan(prod(resolution),2);
                for m = 1:prod(resolution)
                    sol = ode45(@(t,y)derivative(t,y,false),timespan,initialPosition(m,:),odeSolverOptions);
                    finalPositionMainGrid(m,:) = deval(sol,timespan(end));
                end
            end
            
            cgStrainMainGrid = compute_cgStrain(finalPositionMainGrid,initialPosition,resolution);
            
            [~,cgStrainD] = arrayfun(@(x11,x12,x22)eig_array(x11,x12,x22,customEigMethod),cgStrainMainGrid(:,1),cgStrainMainGrid(:,2),cgStrainMainGrid(:,3),'UniformOutput',false);
        end
        
        cgStrainD = cell2mat(cgStrainD);
        
    case 'equationOfVariation'
        
        dFlowMap0 = eye(2);
        dFlowMap0 = reshape(dFlowMap0,4,1);
        nPosition = size(initialPosition,1);
        dFlowMap = nan(nPosition,4);
        
        if coupledIntegration
            initialPosition = [initialPosition,repmat(transpose(dFlowMap0),size(initialPosition,1),1)];
            sol = ode45_vector(@(t,y)derivative(t,y,true),timespan,initialPosition,true,odeSolverOptions);
            dFlowMap = sol(:,3:6);
            % FIXME Check indices in flow definition file.
            dFlowMap(:,[2,3]) = fliplr(dFlowMap(:,[2,3]));
        else
            parfor iPosition = 1:nPosition
                position0 = transpose(initialPosition(iPosition,:));
                y0 = [position0;dFlowMap0];
                sol = ode45(@(t,y)derivative(t,y,true),timespan,y0,odeSolverOptions); %#ok<PFBNS>
                dFlowMap(iPosition,:) = transpose(deval(sol,timespan(end),3:6));
            end
        end
        
        cgStrain = cgStrain_from_dFlowMap(dFlowMap);
        [cgStrainV,cgStrainD] = arrayfun(@(a1,a2,a4)eig_array(a1,a2,a4),squeeze(cgStrain(1,1,:)),squeeze(cgStrain(1,2,:)),squeeze(cgStrain(2,2,:)),'UniformOutput',false);
        cgStrainV = cell2mat(cgStrainV);
        cgStrainD = cell2mat(cgStrainD);
end

if incompressible
    idx = cgStrainD(:,2) < 1;
    n = sum(idx);
    if n
        if n > 1
            warning([mfilename,':incompressible'],['lambda2 < 1 at ',num2str(n),' points; min(lambda2) = ',num2str(min(cgStrainD(idx,2))),'. Eigenvalues and eigenvectors set to NaN at those points.'])
        else
            warning([mfilename,':incompressible'],['lambda2 < 1 at ',num2str(n),' point; lambda2 = ',num2str(cgStrainD(idx,2)),'. Eigenvalues and eigenvectors set to NaN at that point.'])
        end
    end
    prodCgStrainD = prod(cgStrainD(~idx,:),2);
    if any(prodCgStrainD ~= 1)
        warning([mfilename,':prodCgStrainDNotOne'],'min(prodCgStrainD) = %g, max(prodCgStrainD) = %g',min(prodCgStrainD),max(prodCgStrainD))
    end
    cgStrainD(~idx,1) = 1./cgStrainD(~idx,2);
    cgStrainD(idx,:) = nan;
    if nargout == 2
        cgStrainV(idx,:) = nan;
    end
end

switch nargout
    case 1
        cgStrainD = negative_to_nan(cgStrainD);
    case 2
        [cgStrainD,cgStrainV] = negative_to_nan(cgStrainD,cgStrainV);
    otherwise
        error('Number of output arguments incorrect.')
end

switch nargout
    case 1
        varargout{1} = cgStrainD;
    case 2
        varargout{1} = cgStrainV;
        varargout{2} = cgStrainD;
    otherwise
        error('Number of output arguments incorrect.')
end

% eig_array eig function for use with arrayfun
%
% SYNTAX
% [v,d] = eig_array(x11,x12,x22,customMethod)

function [v,d] = eig_array(x11,x12,x22,customMethod)

if customMethod
    [v,d] = eig_custom([x11,x12;x12,x22]);
else    
    [v,d] = eig([x11,x12;x12,x22]);
end

d = transpose(diag(d));
v = reshape(v,1,4);

function cgStrain = cgStrain_from_dFlowMap(dFlowMap)

nRows = size(dFlowMap,1);
dFlowMap = reshape(transpose(dFlowMap),[2 2 nRows]);
cgStrain = arrayfun(@(idx)transpose(dFlowMap(:,:,idx))*dFlowMap(:,:,idx),1:nRows,'UniformOutput',false);
cgStrain = cell2mat(cgStrain);
cgStrain = reshape(cgStrain,[2 2 nRows]);

function [v,d] = eig_custom(a)

d(2,2) = .5*trace(a) + sqrt(.25*trace(a)^2 - det(a));
d(1,1) = .5*trace(a) - sqrt(.25*trace(a)^2 - det(a));
if any(imag(d(:)))
    warning([mfilename,':complexEigenvalue'],['Complex eigenvalue: ',num2str(d([1,4]))])
end

v(1,2) = -a(1,2)/sqrt(a(1,2)^2 + (a(1,1) - d(2,2))^2);
v(2,2) = (a(1,1) - d(2,2))/sqrt(a(1,2)^2 + (a(1,1) - d(2,2))^2);

v(1,1) = v(2,2);
v(2,1) = -v(1,2);

function yf = ode45_vector(odefun,tspan,y0,useEoV,options)

% Reshape m-by-2 array to column array
y0 = transpose(y0);
y0 = y0(:);

if useEoV
    coupledSize = 6;
else
    coupledSize = 2;
end

% Specify three timesteps for ode45's tspan. This has been reported to
% reduce memory usage.
tsteps = [tspan(1),tspan(1) + .5*diff(tspan),tspan(2)];
[~,yf] = ode45(odefun,tsteps,y0,options);
yf = yf(end,:);
yf = transpose(reshape(yf,coupledSize,size(yf,2)/coupledSize));

% Set points with negative eigenvalues to NaN
% The Cauchy-Green strain tensor is positive definite, but numerical
% integration does not enforce this. This function sets those points where
% the Cauchy-Green strain tensor is not positive definite to NaN
function [cgStrainD,varargout] = negative_to_nan(cgStrainD,varargin)

negIdx = any(cgStrainD <= 0,2);

if negIdx
    warning([mfilename,':negativeEigenvalue'],'Negative eigenvalues')
end

cgStrainD(negIdx,:) = nan;

if nargout == 2
    cgStrainV = varargin{1};
    cgStrainV(negIdx,:) = nan;
    varargout{1} = cgStrainV;
end

function auxiliaryPosition = auxiliary_position(basePosition,delta)
%AUXILIARY_POSITION Auxiliary grid positions
%   auxiliaryPosition = AUXILIARY_POSITION(position,delta)
%
%   BASE_POSITION gives coordinates in rows. For 2 dimensional problems
%   it must have two columns and one row for each coordinate. For 3
%   dimensional problems it must have three columns.
%
%   AUXILIARY_POSITION_ has the same format as INITIAL_POSITION, but each
%   row gives positions at auxiliary positions. In 2 dimensions, each
%   row has 8 columns with the format:
%      X+DELTA Y X-DELTA Y X Y+DELTA X Y-DELTA

if size(basePosition,2) == 2
    nColAuxiliaryPosition = 8;
elseif size(basePosition,2) == 3
    nColAuxiliaryPosition = 12;
else
    error('Auxiliary position:Column Number Error')
end

auxiliaryPosition = nan(size(basePosition,1),nColAuxiliaryPosition);

auxiliaryPosition(:,1) = basePosition(:,1) + delta;
auxiliaryPosition(:,2) = basePosition(:,2);
auxiliaryPosition(:,3) = basePosition(:,1) - delta;
auxiliaryPosition(:,4) = basePosition(:,2);
auxiliaryPosition(:,5) = basePosition(:,1);
auxiliaryPosition(:,6) = basePosition(:,2) + delta;
auxiliaryPosition(:,7) = basePosition(:,1);
auxiliaryPosition(:,8) = basePosition(:,2) - delta;

function cgStrain = compute_cgStrain(finalPosition,initialPosition,resolution)
%compute_cgStrain   Compute Cauchy-Green strain
%
%   finalPosition gives positions grouped in rows. In 2 dimensions, each
%   row has 8 columns with the format:
%      (x+delta, y) (x-delta, y) (x, y+delta) (x, y-delta)
%
%   cgStrain gives the Cauchy-Green strain tensor. In 2 dimensions it has
%   the format: c_11 c_12 c_22
% On main grid:
% finalPosition = [x1 y1; x2 y2; ... xN yN]
%
% On auxiliary grid:
% finalPosition = [x1+deltaX y1 x1-deltaX y1 x1 y1+deltaY x1 y1-deltaY;
%                  x2+deltaX y2 x2-deltaX y2 x2 y2+deltaY x2 y2-deltaY;
%                                          ...
%                  xN+deltaX yN xN-deltaX yN xN yN+deltaY xN yN-deltaY]

cgStrain = nan(size(finalPosition,1),3);

switch size(finalPosition,2)
    case 2 % Main grid
        finalX = reshape(finalPosition(:,1),fliplr(resolution));
        finalY = reshape(finalPosition(:,2),fliplr(resolution));
        
        initialPositionX = reshape(initialPosition(:,1),fliplr(resolution));
        initialPositionY = reshape(initialPosition(:,2),fliplr(resolution));
        
        initialPositionX = initialPositionX(1,:);
        initialPositionY = initialPositionY(:,1);
        
        [gradF11,gradF12] = gradient(finalX,initialPositionX,initialPositionY);
        [gradF21,gradF22] = gradient(finalY,initialPositionX,initialPositionY);
        
        gradF11 = reshape(gradF11,prod(double(resolution)),1);
        gradF12 = reshape(gradF12,prod(double(resolution)),1);
        gradF21 = reshape(gradF21,prod(double(resolution)),1);
        gradF22 = reshape(gradF22,prod(double(resolution)),1);
    case 8 % Auxiliary grid
        finalX = finalPosition(:,1:2:7);
        finalY = finalPosition(:,2:2:8);
        
        deltaX = diff(initialPosition(1,[3,1]));
        deltaY = diff(initialPosition(1,[8,6]));
        if deltaX ~= deltaY
            warning([mfilename,':unequalAuxGridDelta'],'Unequal auxiliary grid spacing: (deltaX - deltaY)/min([deltaX,deltaY]) = %.3g.',(deltaX - deltaY)/min([deltaX,deltaY]))
        end
        
        gradF11 = (finalX(:,1) - finalX(:,2))/(deltaX);
        gradF12 = (finalX(:,3) - finalX(:,4))/(deltaY);
        gradF21 = (finalY(:,1) - finalY(:,2))/(deltaX);
        gradF22 = (finalY(:,3) - finalY(:,4))/(deltaY);
        
    otherwise
        error('Number of columns in finalPosition incorrect.')
end

% cgStrain = [c11 c12 c22]
cgStrain(:,1) = gradF11.^2 + gradF21.^2;
cgStrain(:,2) = gradF11.*gradF12 + gradF21.*gradF22;
cgStrain(:,3) = gradF12.^2 + gradF22.^2;
