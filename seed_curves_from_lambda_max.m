% seed_strainlines_from_lambda Seed strainlines from lambda maxima
%
% SYNTAX
% [strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda(distance,cgEigenvalue,cgEigenvector,flowDomain)
% [strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda(distance,cgEigenvalue,cgEigenvector,flowDomain,nMaxStrainlines)
%
% INPUT ARGUMENTS
% distance: threshold distance for placement of lambda_2 maxima
% nMaxStrainlines: Maximum number of hyperbolic LCSs to generate. Default
% is numel(cgEigenvalue).

function [strainlinePosition,strainlineInitialPosition] = seed_strainlines_from_lambda_max(distance,strainlineMaxLength,cgEigenvalue,cgEigenvector,flowDomain,varargin)

narginchk(5,6)

p = inputParser;
addRequired(p,'distance',@(distance)validateattributes(distance,{'double'},{'scalar','>',0}))
addRequired(p,'strainlineMaxLength',@(strainlineMaxLength)validateattributes(strainlineMaxLength,{'double'},{'scalar','>',0}))
addRequired(p,'cgEigenvalue',@(cgEigenvalue)validateattributes(cgEigenvalue,{'double'},{'2d'}))

% Must parse cgEigenvalue before proceeding
parse(p,distance,strainlineMaxLength,cgEigenvalue)
distance = p.Results.distance;
strainlineMaxLength = p.Results.strainlineMaxLength;
cgEigenvalue = p.Results.cgEigenvalue;
flowResolution = fliplr(size(cgEigenvalue));

p = inputParser;
addRequired(p,'cgEigenvector',@(cgEigenvector)validateattributes(cgEigenvector,{'double'},{'size',[fliplr(flowResolution),2]}))
addRequired(p,'flowDomain',@(flowDomain)validateattributes(flowDomain,{'double'},{'size',[2,2]}))
uint = {'uint8','uint16','uint32','uint64'};
addOptional(p,'nMaxStrainlines',numel(cgEigenvalue),@(nMaxStrainlines)validateattributes(nMaxStrainlines,uint,{'scalar','>',0}));

parse(p,cgEigenvector,flowDomain,varargin{:})

cgEigenvector = p.Results.cgEigenvector;
flowDomain = p.Results.flowDomain;
nMaxStrainlines = p.Results.nMaxStrainlines;

%% Compute hyperbolic LCSs seeded from λ₂ local maxima
% Array that records grid points where a strainline already exits
flagArray = false(fliplr(flowResolution));

gridSpace = diff(flowDomain(1,:))/(double(flowResolution(1))-1);
if gridSpace ~= diff(flowDomain(2,:))/(double(flowResolution(2))-1)
    warning([mfilename,':unequalDelta'],['Unequal deltaX (',num2str(gridSpace),') and deltaY (',num2str(diff(flowDomain(2,:))/(double(flowResolution(2))-1)),').'])
end
distanceGridPoints = uint64(distance./gridSpace);

gridPosition{1} = linspace(flowDomain(1,1),flowDomain(1,2),flowResolution(1));
gridPosition{2} = linspace(flowDomain(2,1),flowDomain(2,2),flowResolution(2));

strainlinePosition = cell(1,nMaxStrainlines);
strainlineInitialPosition = nan(2,nMaxStrainlines);

% Find all local maxima with a distance threshold
[cgEigenvalueLocalMax,cgEigenvalueLocalMaxPosition] = local_max2D_gridded(cgEigenvalue,distanceGridPoints);
nMaxStrainlines = min([nMaxStrainlines,numel(cgEigenvalueLocalMax)]);
[cgEigenvalueLocalMax,sortIndex] = sort(cgEigenvalueLocalMax,'descend');
cgEigenvalueLocalMaxPosition = cgEigenvalueLocalMaxPosition(sortIndex,:);

cgEigenvectorInterpolant{1} = griddedInterpolant(fliplr(gridPosition),cgEigenvector(:,:,1));
cgEigenvectorInterpolant{2} = griddedInterpolant(fliplr(gridPosition),cgEigenvector(:,:,2));

% FIXME Should store cgEigenvector as m-by-n array or column array, not
% both
cgEigenvectorColumn = reshape(cgEigenvector,[numel(cgEigenvalue),2]);
nStrainlines = 0;
odeSolverOptions = odeset('relTol',1e-4);

while nStrainlines < nMaxStrainlines
    [nextLocalMax,loc] = find_next_local_max(cgEigenvalueLocalMax,cgEigenvalueLocalMaxPosition,flagArray);
    if isempty(nextLocalMax)
        break
    end
    nStrainlines = nStrainlines + 1;
    strainlineInitialPosition(:,nStrainlines) = [gridPosition{1}(loc(2)),gridPosition{2}(loc(1))];
    periodicBc = false(2,1);
    % Event detection within integrate_line appears unreliable when initial
    % position is on domain boundary, therefore check if
    % strainlineInitialPosition is on boundary before calling
    % integrate_line.
    forwardTime = true;
    backwardTime = true;
    if strainlineInitialPosition(1,nStrainlines) <= flowDomain(1,1)
        if cgEigenvectorInterpolant{1}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) > 0
            backwardTime = false;
        end
        if cgEigenvectorInterpolant{1}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) < 0
            forwardTime = false;
        end
    end
    if strainlineInitialPosition(1,nStrainlines) >= flowDomain(1,2)
        if cgEigenvectorInterpolant{1}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) > 0
            forwardTime = false;
        end
        if cgEigenvectorInterpolant{1}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) < 0
            backwardTime = false;
        end
    end
    if strainlineInitialPosition(2,nStrainlines) <= flowDomain(2,1)
        if cgEigenvectorInterpolant{2}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) < 0
            forwardTime = false;
        end
        if cgEigenvectorInterpolant{2}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) > 0
            backwardTime = false;
        end
    end
    if strainlineInitialPosition(2,nStrainlines) >= flowDomain(2,2)
        if cgEigenvectorInterpolant{2}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) > 0
            forwardTime = false;
        end
        if cgEigenvectorInterpolant{2}(fliplr(transpose(strainlineInitialPosition(:,nStrainlines)))) < 0
            backwardTime = false;
        end
    end
    if forwardTime
        positionPos = integrate_line([0,strainlineMaxLength],strainlineInitialPosition(:,nStrainlines),flowDomain,flowResolution,periodicBc,cgEigenvectorColumn,odeSolverOptions);
    else
        positionPos = [];
    end
    if backwardTime
        positionNeg = integrate_line([0,-strainlineMaxLength],strainlineInitialPosition(:,nStrainlines),flowDomain,flowResolution,periodicBc,cgEigenvectorColumn,odeSolverOptions);
    else
        positionNeg = [];
    end
    % Remove duplicate point from forward and backward time integrations
    if isempty(positionNeg)
        strainlinePosition{nStrainlines} = positionPos;
    else
        strainlinePosition{nStrainlines} = [flipud(positionNeg);positionPos(2:end,:)];
    end
    iFlagArray = line_grid_intersection(strainlinePosition{nStrainlines},gridPosition,gridSpace,distanceGridPoints);
    flagArray = flagArray | iFlagArray;
end

% Remove unused elements
strainlinePosition = strainlinePosition(~cellfun(@(input)isempty(input),strainlinePosition));
strainlineInitialPosition = strainlineInitialPosition(:,1:nStrainlines);

function [nextLocalMax,nextPosition] = find_next_local_max(localMax,position,flagArray)

% Check if next local maximum overlaps with flagArray
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

function flagArray = line_grid_intersection(position,gridPosition,gridSpace,distanceGridPoints)

nX = numel(gridPosition{1});
nY = numel(gridPosition{2});
flagArray = false(nY,nX);

nPoints = size(position,1);

circleMask = circle_mask(distanceGridPoints);
sizeArray = [numel(gridPosition{2}),numel(gridPosition{1})];

offsetGridPosition{1} = gridPosition{1} - .5*gridSpace;
offsetGridPosition{2} = gridPosition{2} - .5*gridSpace;

for iPoint = 1:nPoints-1
    % Create sequence of points along line with spacing no larger than
    % grid spacing.
    deltaX = diff(position([iPoint,iPoint+1],1));
    deltaY = diff(position([iPoint,iPoint+1],2));
    segmentLength = hypot(deltaX,deltaY);
    nPointsSegment = ceil(segmentLength/gridSpace)+1;
    positionSegmentX = linspace(position(iPoint,1),position(iPoint+1,1),nPointsSegment);
    positionSegmentY = linspace(position(iPoint,2),position(iPoint+1,2),nPointsSegment);
    positionSegment = transpose([positionSegmentX;positionSegmentY]);
    for iPointSegment = 1:nPointsSegment
        xI = find(positionSegment(iPointSegment,1) < offsetGridPosition{1},1)-1;
        if isempty(xI)
            xI = nX;
        end
        if xI == 0
            xI = 1;
        end
        yI = find(positionSegment(iPointSegment,2) < offsetGridPosition{2},1)-1;
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
