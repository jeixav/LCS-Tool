function flagArray = line_grid_intersection(position,gridPosition,distanceGridPoints)

nX = numel(gridPosition{1});
nY = numel(gridPosition{2});
flagArray = false(nY,nX);

nPoints = size(position,1);

% showPlot = false;

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
    
    % Loop over circular mask array around each point
%    

%     xI = find(position(iPoint,1) < gridPosition{1},1)-1;
%     if isempty(xI)
%         xI = nX;
%     end
%     if xI == 0
%         xI = 1;
%     end
%     
%     yI = find(position(iPoint,2) < gridPosition{2},1)-1;
%     if isempty(yI)
%         yI = nY;
%     end
%     if yI == 0
%         yI = 1;
%     end
%     
%     centre = [yI,xI];
    
%     circleMaskGlobal = centred_mask(centre,sizeArray,circleMask);
%     flagArray = flagArray | circleMaskGlobal;
    
%     
%     flagArray(yI,xI) = true;
    
%     if showPlot
%         % plot(position(iPoint,1),position(iPoint,2),'*')
%         hImagesc = findobj(gca,'type','image');
%         CData = get(hImagesc,'CData');
%         CData(flagArray) = max(CData(:));
%         set(hImagesc,'CData',CData)
%     end
end
