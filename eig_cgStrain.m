% eig_cgStrain Calculate eigenvalues and eigenvectors of Cauchy-Green strain
%
% SYNTAX
% [cgStrainD,cgStrainV,cgStrain,finalPosition,dFlowMap] = eig_cgStrain(flow,method,customEigMethod,coupledIntegration,verbose)
%
% INPUT ARGUMENTS
% - flow: Structure with fields named derivative, domain, timespan and
% resolution. flow.derivative is a function handle that evaluates the flow
% velocity equations. The function must have the form:
% [velX1;velY1;velX2;velY2;...;velXN;velYN] = f(t,[x1;y1;x2;y2;...;xN;yN].
% flow.domain defines the flow domain. It is of the form:
% [xMin,xMax;yMin,yMax]. flow.timespan defines the flow timespan. It is of
% the form: [tStart,tEnd]. flow.resolution defines the Cauchy-Green strain
% resolution. It is of the form: [resolutionX,resolutionY].
%
% - method: Structure with a field named name that should be either
% 'finiteDifference' or 'equationOfVariation'. Default is finiteDifference.
% If method.name is 'finiteDifference', method.auxiliaryGridRelativeDelta
% can be specified. It should be a number between 0 and 0.5. Default is
% 1e-2. If method.name is 'finiteDifference', method.eigenvalueFromMainGrid
% can be set to true or false to control whether eigenvalues of the
% Cauchy-Green strain are calculated from main grid or auxiliary grid
% points. Default is true.
%
% - customEigMethod: true or false. Default is false.
%
% - coupledIntegration: true or false. Default is true.
%
% - verbose: true or false. Default is false.

function [cgStrainD,cgStrainV,cgStrain,finalPosition,dFlowMap] = eig_cgStrain(flow,varargin)

%% Parse inputs
narginchk(1,5)

p = inputParser;
p.StructExpand = false;
addRequired(p,'flow',@(flow)all(isfield(flow,{'derivative','domain','timespan','resolution'})))
addOptional(p,'method',struct('name','finiteDifference'),@isstruct)
addOptional(p,'customEigMethod',false,@islogical)
addOptional(p,'coupledIntegration',true,@islogical)
addOptional(p,'verbose',false,@islogical)

parse(p,flow,varargin{:})

method = p.Results.method;
customEigMethod = p.Results.customEigMethod;
coupledIntegration = p.Results.coupledIntegration;
verbose = p.Results.verbose;

p = inputParser;
p.KeepUnmatched = true;

addParamValue(p,'odeSolverOptions',[],@isstruct);
parse(p,flow)

odeSolverOptions = p.Results.odeSolverOptions;

%% Main code
initialPosition = initialize_ic_grid(flow.resolution,flow.domain);

switch method.name
    case 'finiteDifference'
        
        %% Parse method structure
        p = inputParser;
        p.KeepUnmatched = true;
        validationFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= .5);
        addParamValue(p,'auxiliaryGridRelativeDelta',1e-2,validationFcn)
        addParamValue(p,'eigenvalueFromMainGrid',true,@islogical)
        parse(p,method)
        
        auxiliaryGridRelativeDelta = p.Results.auxiliaryGridRelativeDelta;
        eigenvalueFromMainGrid = p.Results.eigenvalueFromMainGrid;
        
        %% Eigenvectors from auxiliary grid
        deltaX = (flow.domain(1,2) - flow.domain(1,1))/double(flow.resolution(1))*auxiliaryGridRelativeDelta;
        auxiliaryGridAbsoluteDelta = deltaX;
        auxiliaryPosition = auxiliary_position(initialPosition,auxiliaryGridAbsoluteDelta);
        
        % Transform auxiliaryPosition into a two column array
        auxiliaryPositionX = auxiliaryPosition(:,1:2:end-1);
        auxiliaryPositionY = auxiliaryPosition(:,2:2:end);
        auxiliaryPosition = [auxiliaryPositionX(:) auxiliaryPositionY(:)];
        
        if coupledIntegration
            finalPositionAuxGrid = ode45_vector(@(t,y)flow.derivative(t,y,false),flow.timespan,auxiliaryPosition,false,odeSolverOptions);
        else
            error('Uncoupled integration code not programmed')
        end
        
        % Transform finalPosition into an eight column array
        finalPositionAuxGridX = finalPositionAuxGrid(:,1);
        finalPositionAuxGridY = finalPositionAuxGrid(:,2);
        nPoints = prod(double(flow.resolution));
        finalPositionAuxGridX = reshape(finalPositionAuxGridX,nPoints,4);
        finalPositionAuxGridY = reshape(finalPositionAuxGridY,nPoints,4);
        finalPositionAuxGrid = nan(nPoints,8);
        finalPositionAuxGrid(:,1:2:7) = finalPositionAuxGridX;
        finalPositionAuxGrid(:,2:2:8) = finalPositionAuxGridY;
        
        cgStrainAuxGrid = compute_cgStrain(finalPositionAuxGrid,flow,auxiliaryGridRelativeDelta);
        
        [cgStrainV,cgStrainD] = arrayfun(@(x11,x12,x22)eig_array(x11,x12,x22,customEigMethod),cgStrainAuxGrid(:,1),cgStrainAuxGrid(:,2),cgStrainAuxGrid(:,3),'UniformOutput',false);
        
        cgStrainV = cell2mat(cgStrainV);
        
        %% Eigenvalues from main grid
        if eigenvalueFromMainGrid
            initialPosition = initialize_ic_grid(flow.resolution,flow.domain);
          
            if coupledIntegration
                finalPositionMainGrid = ode45_vector(@(t,y)flow.derivative(t,y,false),flow.timespan,initialPosition,false,odeSolverOptions);
            else
                error('Uncoupled integration code not programmed')
            end
            
            cgStrainMainGrid = compute_cgStrain(finalPositionMainGrid,flow);
            
            [~,cgStrainD] = arrayfun(@eig_array,cgStrainMainGrid(:,1),cgStrainMainGrid(:,2),cgStrainMainGrid(:,3),'UniformOutput',false);
        end
        
        cgStrainD = cell2mat(cgStrainD);
        
        % Use the Cauchy-Green strain calculated with the auxiliary
        % grid for statistics.
        nRows = size(cgStrainAuxGrid,1);
        cgStrain = arrayfun(@(idx)[cgStrainAuxGrid(idx,1),cgStrainAuxGrid(idx,2);cgStrainAuxGrid(idx,2),cgStrainAuxGrid(idx,3)],1:nRows,'uniformOutput',false);
        cgStrain = cell2mat(cgStrain);
        cgStrain = reshape(cgStrain,[2 2 nRows]);
        
    case 'equationOfVariation'
        
        dFlowMap0 = eye(2);
        dFlowMap0 = reshape(dFlowMap0,4,1);
        nPosition = size(initialPosition,1);
        finalPosition = nan(nPosition,2);
        dFlowMap = nan(nPosition,4);
        
        if coupledIntegration
            initialPosition = [initialPosition,repmat(transpose(dFlowMap0),size(initialPosition,1),1)];
            sol = ode45_vector(@(t,y)flow.derivative(t,y,true),flow.timespan,initialPosition,true,odeSolverOptions);
            dFlowMap = sol(:,3:6);
            % FIXME Check indices in flow definition file.
            dFlowMap(:,[2,3]) = fliplr(dFlowMap(:,[2,3]));
            finalPosition = sol(:,1:2);
        else
            parfor iPosition = 1:nPosition
                position0 = transpose(initialPosition(iPosition,:));
                y0 = [position0;dFlowMap0];
                sol = ode45(@(t,y)flow.derivative(t,y,true),flow.timespan,y0,odeSolverOptions); %#ok<PFBNS>
                finalPosition(iPosition,:) = transpose(deval(sol,flow.timespan(end),1:2));
                dFlowMap(iPosition,:) = transpose(deval(sol,flow.timespan(end),3:6));
            end
        end
        
        cgStrain = cgStrain_from_dFlowMap(dFlowMap);
        [cgStrainV,cgStrainD] = arrayfun(@(a1,a2,a4)eig_array(a1,a2,a4),squeeze(cgStrain(1,1,:)),squeeze(cgStrain(1,2,:)),squeeze(cgStrain(2,2,:)),'UniformOutput',false);
        cgStrainV = cell2mat(cgStrainV);
        cgStrainD = cell2mat(cgStrainD);
end

if isfield(flow,'imposeIncompressibility') && flow.imposeIncompressibility == true
    idx = cgStrainD(:,2) < 1;
    n = sum(idx);
    if n
        if n > 1
            warning([mfilename,':imposeIncompressibility'],['Larger eigenvalue less than one at ',num2str(n),' points. Eigenvalues and eigenvectors set to NaN at those points.'])
        else
            warning([mfilename,':imposeIncompressibility'],['Larger eigenvalue less than one at ',num2str(n),' point. Eigenvalues and eigenvectors set to NaN at that point.'])
        end
    end    
    cgStrainD(~idx,1) = 1./cgStrainD(~idx,2);
    cgStrainD(idx,:) = nan;
    cgStrainV(idx,:) = nan;
end

[cgStrainD,cgStrainV,cgStrain] = negative_to_nan(cgStrainD,cgStrainV,cgStrain);

if verbose
    disp('cgStrain_stats:')
    cgStrain_stats(cgStrain,cgStrainV,cgStrainD)
end

% eig_array eig function for use with arrayfun
%
% SYNTAX
% [v,d] = eig_array(x11,x12,x22)
% [v,d] = eig_array(x11,x12,x22,customMethod)

function [v,d] = eig_array(x11,x12,x22,varargin)

narginchk(3,4)

p = inputParser;
addOptional(p,'customMethod',false,@islogical)
parse(p,varargin{:})
customMethod = p.Results.customMethod;

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

sol = ode45(odefun,tspan,y0,options);
yf = deval(sol,tspan(end));
yf = transpose(reshape(yf,coupledSize,size(yf,1)/coupledSize));

% Set points with negative eigenvalues to NaN
% The Cauchy-Green strain tensor is positive definite, but numerical
% integration does not enforce this. This function sets those points where
% the Cauchy-Green strain tensor is not positive definite to NaN
function [cgStrainD,cgStrainV,cgStrain] = negative_to_nan(cgStrainD,cgStrainV,cgStrain)

negIdx = any(cgStrainD <= 0,2);

cgStrainD(negIdx,:) = nan;
cgStrainV(negIdx,:) = nan;
cgStrain(:,:,negIdx) = nan;

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

function cgStrain = compute_cgStrain(finalPosition,flow,auxiliaryGridRelativeDelta)
%compute_cgStrain   Compute Cauchy-Green strain
%   cgStrain = compute_cgStrain(finalPosition,delta)
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
        finalX = reshape(finalPosition(:,1),fliplr(flow.resolution));
        finalY = reshape(finalPosition(:,2),fliplr(flow.resolution));
        
        deltaX = (flow.domain(1,2) - flow.domain(1,1))/double(flow.resolution(1));
        deltaY = (flow.domain(2,2) - flow.domain(2,1))/double(flow.resolution(2));
        
        [gradF11,gradF12] = gradient(finalX,deltaX,deltaY);
        [gradF21,gradF22] = gradient(finalY,deltaX,deltaY);
        
        gradF11 = reshape(gradF11,prod(double(flow.resolution)),1);
        gradF12 = reshape(gradF12,prod(double(flow.resolution)),1);
        gradF21 = reshape(gradF21,prod(double(flow.resolution)),1);
        gradF22 = reshape(gradF22,prod(double(flow.resolution)),1);
    case 8 % Auxiliary grid
        finalX = finalPosition(:,1:2:7);
        finalY = finalPosition(:,2:2:8);
        
        deltaX = diff(flow.domain(1,:))/double(flow.resolution(1)-1)*auxiliaryGridRelativeDelta;
        deltaY = diff(flow.domain(2,:))/double(flow.resolution(2)-1)*auxiliaryGridRelativeDelta;
        if deltaX ~= deltaY
            warning([mfilename,':unequalDelta'],['Unequal deltaX (',num2str(deltaX),') and deltaY (',num2str(deltaY),').'])
        end
        
        gradF11 = (finalX(:,1) - finalX(:,2))/(2*deltaX);
        gradF12 = (finalX(:,3) - finalX(:,4))/(2*deltaY);
        gradF21 = (finalY(:,1) - finalY(:,2))/(2*deltaX);
        gradF22 = (finalY(:,3) - finalY(:,4))/(2*deltaY);
    otherwise
        error('Number of columns in finalPosition incorrect.')
end

% cgStrain = [c11 c12 c22]
cgStrain(:,1) = gradF11.^2 + gradF21.^2;
cgStrain(:,2) = gradF11.*gradF12 + gradF21.*gradF22;
cgStrain(:,3) = gradF12.^2 + gradF22.^2;
