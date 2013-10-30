% seed_curves_from_lambda_max Seed curves from lambda maxima
%
% SYNTAX
% [curvePosition,curveInitialPosition] = seed_curves_from_lambda_max(distance,cgEigenvalue,cgEigenvector,flowDomain)
% [curvePosition,curveInitialPosition] = seed_curves_from_lambda_max(...,'periodicBc',periodicBc)
% [curvePosition,curveInitialPosition] = seed_curves_from_lambda_max(...,'nMaxCurves',nMaxCurves)
% [curvePosition,curveInitialPosition] = seed_curves_from_lambda_max(...,'odeSolverOptions',odeSolverOptions)
%
% INPUT ARGUMENTS
% distance: threshold distance for placement of lambda maxima
% periodicBc: 2-by-1 logical array specifying flow periodic boundary
% conditions. Default is [false,false].
% nMaxCurves: Maximum number of curves to generate. Default is
% numel(cgEigenvalue).
% odeSolverOptions: integrate_line odeSolverOptions input argument

function [curvePosition,curveInitialPosition] = seed_curves_from_lambda_max(distance,curveMaxLength,cgEigenvalue,cgEigenvector,flowDomain,varargin)

if verLessThan('matlab',' 8.1.0')
    % Testing with R2011b shows this function completes without giving
    % an error, but the results are wrong. To avoid such "silent failures",
    % require MATLAB R2013a
    error('MATLAB 8.1.0 or higher is required.')
end

p = inputParser;
addRequired(p,'distance',@(distance)validateattributes(distance,{'double'},{'scalar','>',0}))
addRequired(p,'curveMaxLength',@(curveMaxLength)validateattributes(curveMaxLength,{'double'},{'scalar','>',0}))
addRequired(p,'cgEigenvalue',@(cgEigenvalue)validateattributes(cgEigenvalue,{'double'},{'2d'}))

% Must parse cgEigenvalue before proceeding
parse(p,distance,curveMaxLength,cgEigenvalue)

flowResolution = fliplr(size(cgEigenvalue));

p = inputParser;
addRequired(p,'cgEigenvector',@(cgEigenvector)validateattributes(cgEigenvector,{'double'},{'size',[fliplr(flowResolution),2]}))
addRequired(p,'flowDomain',@(flowDomain)validateattributes(flowDomain,{'double'},{'size',[2,2]}))
addParamValue(p,'periodicBc',[false,false],@(periodicBc)validateattributes(periodicBc,{'logical'},{'size',[1,2]}));
uint = {'uint8','uint16','uint32','uint64'};
addParamValue(p,'nMaxCurves',numel(cgEigenvalue),@(nMaxCurves)validateattributes(nMaxCurves,uint,{'scalar','>',0}));
addParamValue(p,'odeSolverOptions',odeset)

parse(p,cgEigenvector,flowDomain,varargin{:})

periodicBc = p.Results.periodicBc;
nMaxCurves = p.Results.nMaxCurves;
odeSolverOptions = p.Results.odeSolverOptions;

%% Compute hyperbolic LCSs seeded from λ local maxima
% Array that records grid points where a curve already exits
flagArray = false(fliplr(flowResolution));

gridSpace = diff(flowDomain(1,:))/(double(flowResolution(1))-1);
if gridSpace ~= diff(flowDomain(2,:))/(double(flowResolution(2))-1)
    warning([mfilename,':unequalDelta'],['Unequal deltaX (',num2str(gridSpace),') and deltaY (',num2str(diff(flowDomain(2,:))/(double(flowResolution(2))-1)),').'])
end
distanceGridPoints = uint64(distance./gridSpace);

gridPosition{1} = linspace(flowDomain(1,1),flowDomain(1,2),flowResolution(1));
gridPosition{2} = linspace(flowDomain(2,1),flowDomain(2,2),flowResolution(2));

curvePosition = cell(1,nMaxCurves);
curveInitialPosition = nan(2,nMaxCurves);

% Find all local maxima with a distance threshold
[cgEigenvalueLocalMax,cgEigenvalueLocalMaxPosition] = local_max2D_gridded(cgEigenvalue,distanceGridPoints);
nMaxCurves = min([nMaxCurves,numel(cgEigenvalueLocalMax)]);
[cgEigenvalueLocalMax,sortIndex] = sort(cgEigenvalueLocalMax,'descend');
cgEigenvalueLocalMaxPosition = cgEigenvalueLocalMaxPosition(sortIndex,:);

cgEigenvectorInterpolant{1} = griddedInterpolant(fliplr(gridPosition),cgEigenvector(:,:,1));
cgEigenvectorInterpolant{2} = griddedInterpolant(fliplr(gridPosition),cgEigenvector(:,:,2));

% FIXME Should store cgEigenvector as m-by-n array or column array, not
% both
cgEigenvectorColumn = reshape(cgEigenvector,[numel(cgEigenvalue),2]);
nCurves = 0;

while nCurves < nMaxCurves
    [nextLocalMax,loc] = find_next_local_max(cgEigenvalueLocalMax,cgEigenvalueLocalMaxPosition,flagArray);
    if isempty(nextLocalMax)
        break
    end
    nCurves = nCurves + 1;
    curveInitialPosition(:,nCurves) = [gridPosition{1}(loc(2)),gridPosition{2}(loc(1))];
    % Event detection within integrate_line appears unreliable when initial
    % position is on domain boundary, therefore check if
    % curveInitialPosition is on boundary before calling integrate_line.
    forwardTime = true;
    backwardTime = true;
    if curveInitialPosition(1,nCurves) <= flowDomain(1,1)
        if cgEigenvectorInterpolant{1}(fliplr(transpose(curveInitialPosition(:,nCurves)))) > 0
            backwardTime = false;
        end
        if cgEigenvectorInterpolant{1}(fliplr(transpose(curveInitialPosition(:,nCurves)))) < 0
            forwardTime = false;
        end
    end
    if curveInitialPosition(1,nCurves) >= flowDomain(1,2)
        if cgEigenvectorInterpolant{1}(fliplr(transpose(curveInitialPosition(:,nCurves)))) > 0
            forwardTime = false;
        end
        if cgEigenvectorInterpolant{1}(fliplr(transpose(curveInitialPosition(:,nCurves)))) < 0
            backwardTime = false;
        end
    end
    if curveInitialPosition(2,nCurves) <= flowDomain(2,1)
        if cgEigenvectorInterpolant{2}(fliplr(transpose(curveInitialPosition(:,nCurves)))) < 0
            forwardTime = false;
        end
        if cgEigenvectorInterpolant{2}(fliplr(transpose(curveInitialPosition(:,nCurves)))) > 0
            backwardTime = false;
        end
    end
    if curveInitialPosition(2,nCurves) >= flowDomain(2,2)
        if cgEigenvectorInterpolant{2}(fliplr(transpose(curveInitialPosition(:,nCurves)))) > 0
            forwardTime = false;
        end
        if cgEigenvectorInterpolant{2}(fliplr(transpose(curveInitialPosition(:,nCurves)))) < 0
            backwardTime = false;
        end
    end
    if forwardTime
        positionPos = integrate_line([0,curveMaxLength],curveInitialPosition(:,nCurves),flowDomain,flowResolution,periodicBc,cgEigenvectorColumn,odeSolverOptions);
    else
        positionPos = [];
    end
    if backwardTime
        positionNeg = integrate_line([0,-curveMaxLength],curveInitialPosition(:,nCurves),flowDomain,flowResolution,periodicBc,cgEigenvectorColumn,odeSolverOptions);
    else
        positionNeg = [];
    end
    % Remove duplicate point from forward and backward time integrations
    if isempty(positionNeg)
        curvePosition{nCurves} = positionPos;
    else
        curvePosition{nCurves} = [flipud(positionNeg);positionPos(2:end,:)];
    end
    iFlagArray = line_grid_intersection(curvePosition{nCurves},gridPosition,gridSpace,distanceGridPoints);
    flagArray = flagArray | iFlagArray;
end

% Remove unused cell array elements
curvePosition = curvePosition(~cellfun(@(input)isempty(input),curvePosition));
curveInitialPosition = curveInitialPosition(:,1:nCurves);

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
        % Discard NaN elements as potential local maximums
        % FIXME Would prefer not having to deal with NaNs in array.
        if ~isnan(array(m,n))
            centredMask = centred_mask([m,n],size(array),circle);
            if ~any(array(m,n) < array(centredMask))
                c(m,n) = array(m,n);
                i(m,n) = true;
            end
        end
    end
end

c = c(~isnan(c));
[row,col] = find(i);
i = [row,col];

% Enlarge a circular mask array to an array of size sizeArray centered at
% centre. This involves "padding" circleArray with zeros.
function centredMask = centred_mask(centre,sizeArray,circleMask)

% validateattributes(sizeArray,uint,{'size',[1,2],'>=',1},mfilename,'sizeArray',2)
% validateattributes(centre,uint,{'size',[1,2],'>=',1},mfilename,'centre',1)
% validateattributes(centre(1),uint,{'scalar','<=',sizeArray(1)},mfilename,'centre(1)',1)
% validateattributes(centre(2),uint,{'scalar','<=',sizeArray(2)},mfilename,'centre(2)',1)
% validateattributes(circleMask,{'logical'},{'square'})

if ~rem(size(circleMask,1),2)
    error('circleMask must have odd dimensions')
else
    radius = (size(circleMask,1) - 1)/2 + 1;
end

% Enlarged array notation
% ┌─┬─┬─┐
% │ │1│ │
% │ ├─┤ │
% │4│K│2│
% │ ├─┤ │
% │ │3│ │
% └─┴─┴─┘

% Enlarge array with false entries
array = true(sizeArray);

enlarge1nRows = radius - centre(1);
if enlarge1nRows < 0
    enlarge1nRows = 0;
end
enlarge1nCols = sizeArray(2);
enlarge1 = false(enlarge1nRows,enlarge1nCols);

enlarge3nRows = centre(1) + radius - 1 - sizeArray(1);
if enlarge3nRows < 0
    enlarge3nRows = 0;
end
enlarge3nCols = sizeArray(2);
enlarge3 = false(enlarge3nRows,enlarge3nCols);

enlarge2nRows = sizeArray(1) + enlarge1nRows + enlarge3nRows;
enlarge2nCols = centre(2) + radius - 1 - sizeArray(2);
enlarge2 = false(enlarge2nRows,enlarge2nCols);

enlarge4nRows = sizeArray(1) + enlarge1nRows + enlarge3nRows;
enlarge4nCols = radius - centre(2);
enlarge4 = false(enlarge4nRows,enlarge4nCols);

enlargedArray = [enlarge4,[enlarge1;array;enlarge3],enlarge2];

% Enlarge circle mask with false entries
enlarge1nRows = centre(1) - radius;
if enlarge1nRows < 0
    enlarge1nRows = 0;
end
enlarge1nCols = size(circleMask,1);
enlarge1 = false(enlarge1nRows,enlarge1nCols);

enlarge3nRows = sizeArray(1) - (radius - 1) - centre(1);
if enlarge3nRows < 0
    enlarge3nRows = 0;
end
enlarge3nCols = size(circleMask,1);
enlarge3 = false(enlarge3nRows,enlarge3nCols);

enlarge2nRows = size(circleMask,1) + enlarge1nRows + enlarge3nRows;
enlarge2nCols = sizeArray(2) - radius - centre(2) + 1;
enlarge2 = false(enlarge2nRows,enlarge2nCols);

enlarge4nRows = size(circleMask,1) + enlarge1nRows + enlarge3nRows;
enlarge4nCols = centre(2) - radius;
enlarge4 = false(enlarge4nRows,enlarge4nCols);

enlargedCircle = [enlarge4,[enlarge1;circleMask;enlarge3],enlarge2];

centredMask = enlargedArray & enlargedCircle;
corner1 = [radius-centre(1)+1,radius-centre(2)+1];
corner1(corner1 < 1) = 1;
corner2 = corner1 + sizeArray - [1,1];
centredMask = centredMask(corner1(1):corner2(1),corner1(2):corner2(2));

function circle = circle_mask(radius)

[X,Y] = meshgrid(0:radius,radius:-1:0);
distance = sqrt(double(X).^2 + double(Y).^2);
quarterCircle = distance <= radius;
circle = [[fliplr(quarterCircle),quarterCircle(:,2:end)];[rot90(quarterCircle(1:end-1,:),2),flipud(quarterCircle(1:end-1,2:end))]];
