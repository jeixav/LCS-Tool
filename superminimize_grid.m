function superminIndex = superminimize_grid(position,segmentIndex,...
    relativeStretching,superminDistance,flowDomain,superminResolution,...
    showPlotSuperminLine,superminMethod)

verbose.progress = true;

if nargin < 7
    showPlotSuperminLine = false;
else
    if ishandle(showPlotSuperminLine)
        tmpAxes = showPlotSuperminLine;
        clear('showPlotSuperminLine')
        showPlotSuperminLine.axes = tmpAxes;
        clear('tmpAxes')
        showPlotSuperminLine.flowDomain = flowDomain;
    end
end

% deltaX = (flowDomain(1,2) - flowDomain(1,1))/(strainlineResolution(1) + 1);
% x = flowDomain(1,1) + deltaX:deltaX:flowDomain(1,2) - deltaX;
tmp = initialize_ic_grid(superminResolution,flowDomain);
tmp = reshape(tmp(:,1),fliplr(superminResolution));
x = tmp(1,:);

if verbose.progress
    progressBar = ParforProgressStarter2(mfilename,...
        sum(superminResolution));
end
superminIndexArray = cell(1,superminResolution(1));
parfor i = 1:superminResolution(1)
    superminIndexArray{i} = superminArray(x(i),'x',position,...
        segmentIndex,relativeStretching,superminDistance,...
        showPlotSuperminLine,superminMethod);
    if verbose.progress %#ok<PFBNS>
        progressBar.increment(i) %#ok<PFBNS>
    end
end
superminIndex = superminIndexArray{1};

for m = 2:size(superminIndexArray,2)
    for n = 1:size(superminIndexArray{m},2)
        superminIndex{n} = any([superminIndex{n} superminIndexArray{m}{n}],2);
    end
end

% deltaY = (flowDomain(2,2) - flowDomain(2,1))/(strainlineResolution(2) + 1);
% y = flowDomain(2,1) + deltaY:deltaY:flowDomain(2,2) - deltaY;
tmp = initialize_ic_grid(superminResolution,flowDomain);
tmp = reshape(tmp(:,2),fliplr(superminResolution));
y = tmp(:,1);

superminIndexArray = cell(1,superminResolution(2));
parfor i = 1:superminResolution(2)
    superminIndexArray{i} = superminArray(y(i),'y',position,...
        segmentIndex,relativeStretching,superminDistance,...
        showPlotSuperminLine,superminMethod);
    if verbose.progress %#ok<PFBNS>
        progressBar.increment(i) %#ok<PFBNS>
    end
end
if verbose.progress
    try
        delete(progressBar);
    catch me %#ok<NASGU>
    end
end

for m = 1:size(superminIndexArray,2)
    for n = 1:size(superminIndexArray{m},2)
        superminIndex{n} = any([superminIndex{n} superminIndexArray{m}{n}],2);
    end
end

end

function superminIndexArray = superminArray(xOrY,xOrYFlag,position,...
    segmentIndex,relativeStretching,superminDistance,...
    showPlotSuperminLine,superminMethod)

if nargin < 7
    showPlotSuperminLine = false;
end

superminIndexArray = cellfun(@(x) ~boolean(x),relativeStretching,...
    'UniformOutput',false);

if strcmp(xOrYFlag,'x')
    linePosition = [xOrY 0; xOrY 1];
else
    linePosition = [0 xOrY; 1 xOrY];
end

intersectionPosition = [];
relativeStretchingLine = [];
indexLine = [];

nStrainlines = size(position,2);

for iStrainline = 1:nStrainlines

    positionLocal = position{iStrainline};
    
    nSegments = size(segmentIndex{iStrainline},1);

    for iSegment = 1:nSegments

        segmentPosition = positionLocal(...
            segmentIndex{iStrainline}(iSegment,1):...
            segmentIndex{iStrainline}(iSegment,2),:);
        
%         if isa(showPlotSuperminLine,'numeric')
%             plot(segmentPosition(:,1),segmentPosition(:,2),'r',...
%                 'Tag','segment','linewidth',2)
%         end
        
        sidePointLine = side_point_line(segmentPosition,linePosition);
        
        % Segment point lies exactly on superminimization line.
        segmentPointOnLineIndex = find(sidePointLine == 0);
        if ~isempty(segmentPointOnLineIndex)
            % FIXME Error if numel(segmentPointOnLineIndex) > 1
            if strcmp(xOrYFlag,'x')
                intersectionPosition(end+1) = ...
                    segmentPosition(segmentPointOnLineIndex,2); %#ok<AGROW>
            else
                intersectionPosition(end+1) = ...
                    segmentPosition(segmentPointOnLineIndex,1); %#ok<AGROW>
            end
            relativeStretchingLine(end+1) = ...
                relativeStretching{iStrainline}(iSegment); %#ok<AGROW>
            indexLine(1:2,end+1) = [iStrainline; iSegment]; %#ok<AGROW>
        end
        
        % Segment points lie on both sides of superminimization line.
        sideIndex = find(abs(diff(sidePointLine)) == 2);
        
        if ~isempty(sideIndex)
            for i = transpose(sideIndex)
                % Generalize to have common code for x and y directions.
                % http://www.mathworks.com/matlabcentral/fileexchange/35606-line-line-intersection-2d
                if strcmp(xOrYFlag,'x')
                    intersectionPosition(end+1) = ...
                        interp1(segmentPosition(i:i+1,1),...
                        segmentPosition(i:i+1,2),xOrY); %#ok<AGROW>
                else
                    intersectionPosition(end+1) = ...
                        interp1(segmentPosition(i:i+1,2),...
                        segmentPosition(i:i+1,1),xOrY); %#ok<AGROW>
                end                    
                relativeStretchingLine(end+1) = ...
                    relativeStretching{iStrainline}(iSegment); %#ok<AGROW>
                indexLine(1:2,end+1) = [iStrainline; iSegment]; %#ok<AGROW>
            end
        end
     end
end

% if isa(showPlotSuperminLine,'struct')
    % superminFigure = figure;
    % if strcmp(xOrYFlag,'x')
    %    xLim = showPlotSuperminLine.flowDomain(2,:);
    % else
    %    xLim = showPlotSuperminLine.flowDomain(1,:);
    % end
    % superminAxes = axes('nextplot','add','box','on',...
    %    'xgrid','on','ygrid','on','xlim',xLim,'parent',superminFigure);
    % text(.9,-.1,0,[xOrYFlag,' = ',num2str(xOrY)],'parent',...
    %    superminAxes,'units','normalized')

    % Highlight superminimization segments
    
    % position{indexLine(1,1)}(segmentIndex{indexLine(1,1)}(indexLine(1,2),:),1)
    % position{indexLine(1,1)}(segmentIndex{indexLine(1,1)}(indexLine(1,2),:),1)
% else
%     showPlotSuperminLine = false;
% end

% if isa(showPlotSuperminLine,'struct')
switch superminMethod
    case 'superminimize'
        superminIndexLine = superminimize(intersectionPosition,...
            relativeStretchingLine,superminDistance,...
            showPlotSuperminLine);
    case 'superminimize2'
        superminIndexLine = superminimize2(intersectionPosition,...
            relativeStretchingLine,superminDistance,...
            showPlotSuperminLine);
    otherwise
        error('superminMethod invalid')
end
% else
%     superminIndexLine = superminimize2(intersectionPosition,...
%         relativeStretchingLine,superminDistance);
% end

nSupermin = length(superminIndexLine);

% Highlight superminimzed segments
if isa(showPlotSuperminLine,'struct')
    iSuperminStrainline = indexLine(1,superminIndexLine);
    iSuperminSegment = indexLine(2,superminIndexLine);
    
    superminStartIndex = arrayfun(@(idx)segmentIndex{iSuperminStrainline(idx)}(iSuperminSegment(idx),1),1:nSupermin);
    superminEndIndex = arrayfun(@(idx)segmentIndex{iSuperminStrainline(idx)}(iSuperminSegment(idx),2),1:nSupermin);
    superminPosition = arrayfun(...
        @(idx)position{iSuperminStrainline(idx)}...
        (superminStartIndex(idx):superminEndIndex(idx),:),1:nSupermin,...
        'UniformOutput',false);
    % arrayfun(@(idx) plotSupermin(showPlotSuperminLine.axes,...
    %     superminPosition{idx}),1:nSupermin,'UniformOutput',false)
    % pause
    % delete(superminFigure)
end

for i = 1:nSupermin
    superminIndexArray{indexLine(1,superminIndexLine(i))}...
        (indexLine(2,superminIndexLine(i))) = true;
end

end

function plotSupermin(axes,superminPosition)

plot(axes,superminPosition(:,1),superminPosition(:,2),...
        'Color','r','linewidth',2,'Tag','superminSegment')
    
end