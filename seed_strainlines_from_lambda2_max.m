% seed_strainlines_from_lambda2_max Seeding strainlines from lambda_2 maxima
%
% SYNTAX
% [strainlinePosition,hFigure] = seed_strainlines_from_lambda2_max(distance)
% [strainlinePosition,hFigure] = seed_strainlines_from_lambda2_max(distance,hyperbolicLcsMaxNo)
%
% INPUT ARGUMENTS
% distance: threshold distance for placement of lambda_2 maxima
% hyperbolicLcsMaxNo: Maximum number of hyperbolic LCSs to generate. If not
% specified, this number defaults to the flow resolution.
%
% EXAMPLES
% strainlinePosition = seed_strainlines_from_lambda2_max(.025);
% strainlinePosition = seed_strainlines_from_lambda2_max(.025,uint8(5));

function [strainlinePosition,hFigure] = seed_strainlines_from_lambda2_max(distance,varargin)

narginchk(1,2)

p = inputParser;
addRequired(p,'distance',@(distance)validateattributes(distance,{'double'},{'scalar','>',0}))
uint = {'uint8','uint16','uint32','uint64'};
doubleGyre = double_gyre;
flow = doubleGyre.flow;
flow = set_flow_resolution([251,126],flow);
flow = set_flow_timespan([0,20],flow);
flow = set_flow_ode_solver_options(odeset('relTol',1e-4),flow);

% domainEnlargement = .05;
% flow = set_flow_domain([[0,2]+diff([0,2])*domainEnlargement*[-1,1];[0,1]+diff([0,1])*domainEnlargement*[-1,1]],flow);

addOptional(p,'hyperbolicLcsMaxNo',flow.resolution(1)*flow.resolution(2),@(hyperbolicLcsMaxNo)validateattributes(hyperbolicLcsMaxNo,uint,{'scalar','>',0}));

parse(p,distance,varargin{:})
distance = p.Results.distance;
hyperbolicLcsMaxNo = p.Results.hyperbolicLcsMaxNo;

strainline = doubleGyre.strainline;
clear('doubleGyre')

%% Compute λ₂ array
method.name = 'finiteDifference';
customEigMethod = false;
nBlock = 1;
coupledIntegration = flow.resolution(1)*flow.resolution(2)*4*2/nBlock;
[cgEigenvalue,cgEigenvector] = eig_cgStrain(flow,method,customEigMethod,coupledIntegration);

%% Plot λ₂ field
lambda2 = reshape(cgEigenvalue(:,2),fliplr(flow.resolution));
hFigure = figure;
hAxes = axes;
set(hAxes,'parent',hFigure)
set(hAxes,'NextPlot','add')
set(hAxes,'box','on')
set(hAxes,'DataAspectRatio',[1,1,1])
hImagesc = imagesc(flow.domain(1,:),flow.domain(2,:),log(lambda2));
set(hImagesc,'parent',hAxes)
axis(hAxes,'tight')
hColorbar = colorbar('peer',hAxes);
set(get(hColorbar,'xlabel'),'string','log(\lambda_2)')
drawnow

% %% Plot ξ₁ field
% xPos = linspace(flow.domain(1,1),flow.domain(1,2),flow.resolution(1));
% yPos = linspace(flow.domain(2,1),flow.domain(2,2),flow.resolution(2));
% quiver(xPos,yPos,reshape(cgEigenvector(:,1),fliplr(flow.resolution)),reshape(cgEigenvector(:,2),fliplr(flow.resolution)))

%% Compute hyperbolic LCSs seeded from λ₂ local maxima

% Array that records grid points where a strainline already exits
flagArray = false(size(lambda2));

deltaX = diff(flow.domain(1,:))/double(flow.resolution(1)-1);
deltaY = diff(flow.domain(2,:))/double(flow.resolution(2)-1);
xPos = linspace(flow.domain(1,1),flow.domain(1,2),flow.resolution(1)) - .5*deltaX;
yPos = linspace(flow.domain(2,1),flow.domain(2,2),flow.resolution(2)) - .5*deltaY;
gridPosition{1} = xPos;
gridPosition{2} = yPos;

cmap = colormap(hAxes);
cmap(end,:) = ones(3,1);
colormap(hAxes,cmap)

if deltaX ~= deltaY
    error('Cannot set distnace in units of grid points if deltaX ~= deltaY')
else
    distanceGridPoints = uint64(distance./deltaX);
end

strainlinePosition = cell(1,hyperbolicLcsMaxNo);

nHyperbolicLcs = 0;
odeSolverOptions = odeset('relTol',1e-4);
while nHyperbolicLcs < hyperbolicLcsMaxNo
    [nextLocalMax,loc] = find_next_local_max(lambda2,flagArray,distanceGridPoints);
    if isempty(nextLocalMax)
        break
    end
    initialPosition = [xPos(loc(2))+.5*deltaX,yPos(loc(1))+.5*deltaX];
    hInitialPosition = plot(initialPosition(1),initialPosition(2));
    set(hInitialPosition,'marker','o')
    set(hInitialPosition,'markerEdgeColor','k')
    set(hInitialPosition,'markerFaceColor','k')
    set(hInitialPosition,'markerSize',8)
    periodicBc = false(2,1);
    positionPos = integrate_line([0,strainline.maxLength],initialPosition,flow.domain,flow.resolution,periodicBc,cgEigenvector(:,1:2),odeSolverOptions);
    positionNeg = integrate_line([0,-strainline.maxLength],initialPosition,flow.domain,flow.resolution,periodicBc,cgEigenvector(:,1:2),odeSolverOptions);
    strainlinePosition{nHyperbolicLcs+1} = [flipud(positionNeg);positionPos(2:end,:)];
    hStrainline = plot(hAxes,strainlinePosition{nHyperbolicLcs+1}(:,1),strainlinePosition{nHyperbolicLcs+1}(:,2));
    set(hStrainline,'color','k')
    set(hStrainline,'lineWidth',2)
        
    % Flag grid cells intersection by strainline
    iFlagArray = line_grid_intersection(strainlinePosition{nHyperbolicLcs+1},gridPosition,distanceGridPoints);
    flagArray = flagArray | iFlagArray;
    
    CData = get(findobj(hAxes,'Type','image'),'CData');
    CData(flagArray) = max(CData(:));
    set(findobj(hAxes,'Type','image'),'CData',CData)
    drawnow
    
    nHyperbolicLcs = nHyperbolicLcs + 1;
end

% Remove empty elements of cell array
strainlinePosition = strainlinePosition(~cellfun(@(input)isempty(input),strainlinePosition));
