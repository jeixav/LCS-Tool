% seed_strainlines_from_lambda2_max Seeding strainlines from lambda_2 maxima
%
% SYNTAX
% [strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda2_max(distance,cgEigenvalue2,cgEigenvector1,flowDomain)
% [strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda2_max(distance,cgEigenvalue2,cgEigenvector1,flowDomain,hyperbolicLcsMaxNo)
%
% INPUT ARGUMENTS
% distance: threshold distance for placement of lambda_2 maxima
% hyperbolicLcsMaxNo: Maximum number of hyperbolic LCSs to generate.
% Default is numel(cgEigenvalue2).

function [strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda2_max(distance,strainlineMaxLength,cgEigenvalue2,cgEigenvector1,flowDomain,varargin)

narginchk(4,5)

p = inputParser;
addRequired(p,'distance',@(distance)validateattributes(distance,{'double'},{'scalar','>',0}))
addRequired(p,'strainlineMaxLength',@(strainlineMaxLength)validateattributes(strainlineMaxLength,{'double'},{'scalar','>',0}))
addRequired(p,'cgEigenvalue2',@(cgEigenvalue2)validateattributes(cgEigenvalue2,{'double'},{'2d'}))

% Must parse cgEigenvalue2 before proceeding
parse(p,distance,strainlineMaxLength,cgEigenvalue2)
distance = p.Results.distance;
strainlineMaxLength = p.Results.strainlineMaxLength;
cgEigenvalue2 = p.Results.cgEigenvalue2;
flowResolution = fliplr(size(cgEigenvalue2));

p = inputParser;
addRequired(p,'cgEigenvector1',@(cgEigenvector1)validateattributes(cgEigenvector1,{'double'},{'size',[fliplr(flowResolution),2]}))
addRequired(p,'flowDomain',@(flowDomain)validateattributes(flowDomain,{'double'},{'size',[2,2]}))
uint = {'uint8','uint16','uint32','uint64'};
addOptional(p,'hyperbolicLcsMaxNo',numel(cgEigenvalue2),@(hyperbolicLcsMaxNo)validateattributes(hyperbolicLcsMaxNo,uint,{'scalar','>',0}));

parse(p,cgEigenvector1,flowDomain,varargin{:})

cgEigenvector1 = p.Results.cgEigenvector1;
flowDomain = p.Results.flowDomain;
hyperbolicLcsMaxNo = p.Results.hyperbolicLcsMaxNo;

%% Compute hyperbolic LCSs seeded from λ₂ local maxima
% Array that records grid points where a strainline already exits
flagArray = false(fliplr(flowResolution));

deltaX = diff(flowDomain(1,:))/double(flowResolution(1)-1);
deltaY = diff(flowDomain(2,:))/double(flowResolution(2)-1);
xPos = linspace(flowDomain(1,1),flowDomain(1,2),flowResolution(1)) - .5*deltaX;
yPos = linspace(flowDomain(2,1),flowDomain(2,2),flowResolution(2)) - .5*deltaY;
gridPosition{1} = xPos;
gridPosition{2} = yPos;

if deltaX ~= deltaY
    error(['Cannot set distance in units of grid points because deltaX ~= deltaY. deltaX = ',num2str(deltaX),' deltaY = ',num2str(deltaY)])
else
    distanceGridPoints = uint64(distance./deltaX);
end

strainlinePosition = cell(1,hyperbolicLcsMaxNo);
strainlineInitialPosition = nan(2,hyperbolicLcsMaxNo);

cgEigenvector1 = reshape(cgEigenvector1,[numel(cgEigenvalue2),2]);
nHyperbolicLcs = 0;
odeSolverOptions = odeset('relTol',1e-4);
while nHyperbolicLcs < hyperbolicLcsMaxNo
    [nextLocalMax,loc] = find_next_local_max(cgEigenvalue2,flagArray,distanceGridPoints);
    if isempty(nextLocalMax)
        break
    end
    strainlineInitialPosition(:,nHyperbolicLcs+1) = [xPos(loc(2))+.5*deltaX,yPos(loc(1))+.5*deltaX];
    periodicBc = false(2,1);
    positionPos = integrate_line([0,strainlineMaxLength],strainlineInitialPosition(:,nHyperbolicLcs+1),flowDomain,flowResolution,periodicBc,cgEigenvector1,odeSolverOptions);
    positionNeg = integrate_line([0,-strainlineMaxLength],strainlineInitialPosition(:,nHyperbolicLcs+1),flowDomain,flowResolution,periodicBc,cgEigenvector1,odeSolverOptions);
    strainlinePosition{nHyperbolicLcs+1} = [flipud(positionNeg);positionPos(2:end,:)];
    iFlagArray = line_grid_intersection(strainlinePosition{nHyperbolicLcs+1},gridPosition,distanceGridPoints);
    flagArray = flagArray | iFlagArray;
    nHyperbolicLcs = nHyperbolicLcs + 1;
end

% Remove unused elements
strainlinePosition = strainlinePosition(~cellfun(@(input)isempty(input),strainlinePosition));
strainlineInitialPosition = strainlineInitialPosition(:,1:nHyperbolicLcs);

function [nextLocalMax,nextPosition] = find_next_local_max(array,flagArray,distance)

%% Find all local maxima with a distance threshold
[localMax,position] = local_max2D_gridded(array,distance);
[localMax,sortIndex] = sort(localMax,'descend');
position = position(sortIndex,:);

%% Check if next local maximum overlaps with flagArray
idx = 1;
while idx <= numel(localMax)
    if ~flagArray(position(idx,1),position(idx,2))
        nextLocalMax = localMax(idx);
        nextPosition = position(idx,:);
        return
    end
    idx = idx + 1;
end

% No local maxima left; return empty arrays
nextLocalMax = [];
nextPosition = [];

function flagArray = line_grid_intersection(position,gridPosition,distanceGridPoints)

nX = numel(gridPosition{1});
nY = numel(gridPosition{2});
flagArray = false(nY,nX);

nPoints = size(position,1);

circleMask = circle_mask(distanceGridPoints);
sizeArray = [numel(gridPosition{2}),numel(gridPosition{1})];

gridDeltaX = mean(diff(gridPosition{1}));
gridDeltaY = mean(diff(gridPosition{2}));
if ~(gridDeltaX == gridDeltaY)
    warning([mfilename,':unequalGridSpace'],['Unequal deltaX (',num2str(gridDeltaX),') and deltaY (',num2str(gridDeltaY),').'])
end

for iPoint = 1:nPoints-1
    % Create sequence of points along line with spacing no larger than
    % grid spacing.
    deltaX = diff(position([iPoint,iPoint+1],1));
    deltaY = diff(position([iPoint,iPoint+1],2));
    segmentLength = hypot(deltaX,deltaY);
    nPointsSegment = ceil(segmentLength/gridDeltaX)+1;
    positionSegmentX = linspace(position(iPoint,1),position(iPoint+1,1),nPointsSegment);
    positionSegmentY = linspace(position(iPoint,2),position(iPoint+1,2),nPointsSegment);
    positionSegment = transpose([positionSegmentX;positionSegmentY]);
    for iPointSegment = 1:nPointsSegment
        xI = find(positionSegment(iPointSegment,1) < gridPosition{1},1)-1;
        if isempty(xI)
            xI = nX;
        end
        if xI == 0
            xI = 1;
        end
        yI = find(positionSegment(iPointSegment,2) < gridPosition{2},1)-1;
        if isempty(yI)
            yI = nY;
        end
        if yI == 0
            yI = 1;
        end
        centre = [yI,xI];
        circleMaskGlobal = centred_mask(centre,sizeArray,circleMask);
        flagArray = flagArray | circleMaskGlobal;
    end
end

% Locally largest element in 2D array.
%
% SYNTAX
% [c,i] = local_max2D_gridded(array,radius)
%
% EXAMPLE
% C = peaks;
% [c,i] = local_max2D_gridded(C,uint8(5));
% imagesc(C)
% hold on
% plot(i(:,2),i(:,1),'ko')

function [c,i] = local_max2D_gridded(array,radius)

validateattributes(array,{'numeric'},{'2d'})
validateattributes(radius,{'uint8','uint16','uint32','uint64'},{'scalar','>',0})

circle = circle_mask(radius);

c = nan(size(array));
i = false(size(array));

for m = 1:size(array,1)
    for n = 1:size(array,2)
        centredMask = centred_mask([m,n],size(array),circle);
        if ~any(array(m,n) < array(centredMask))
            c(m,n) = array(m,n);
            i(m,n) = true;
        end
    end
end

c = c(~isnan(c));
[row,col] = find(i);
i = [row,col];
