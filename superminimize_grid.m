% superminimize_grid Superminimization on a Cartesian grid of lines
%
% EXAMPLE
% To see a graphical representation of the superminimization process, set
% verbose.graphs = true. For example:
%
% matlabpool('open')
% pctRunOnAll javaaddpath('ParforProgress2')
% addpath('flow_templates')
% bickleyJet = bickley_jet(1);
% verbose.graphs = true;
% bickleyJet = strain_lcs_script(bickleyJet,[],verbose);

function superminIndex = superminimize_grid(position,segmentIndex,...
    relativeStretching,superminDistance,flowDomain,superminResolution,...
    showPlotSuperminLine,verbose)

if nargin < 7
    showPlotSuperminLine = false;
end

if showPlotSuperminLine
    hAxes = setup_figure(flowDomain);

    hStrainline = cellfun(@(position)plot(hAxes,position(:,1),...
        position(:,2)),position);
    set(hStrainline,'tag','strainline')
    set(hStrainline,'color',.8*[1 1 1])
    showPlotSuperminLine = hAxes;
end

tmp = initialize_ic_grid(superminResolution,flowDomain);
tmp = reshape(tmp(:,1),fliplr(superminResolution));
x = tmp(1,:);

if verbose
    progressBar = ParforProgressStarter2(mfilename,...
        sum(superminResolution));
else
    progressBar = [];
end
superminIndexArray = cell(1,superminResolution(1));
parfor i = 1:superminResolution(1)
    superminIndexArray{i} = superminArray(x(i),'x',position,...
        segmentIndex,relativeStretching,superminDistance,...
        showPlotSuperminLine);
    if verbose 
        progressBar.increment(i) %#ok<PFBNS>
    end
end
superminIndex = superminIndexArray{1};

for m = 2:size(superminIndexArray,2)
    for n = 1:size(superminIndexArray{m},2)
        superminIndex{n} = any([superminIndex{n} superminIndexArray{m}{n}],2);
    end
end

tmp = initialize_ic_grid(superminResolution,flowDomain);
tmp = reshape(tmp(:,2),fliplr(superminResolution));
y = tmp(:,1);

superminIndexArray = cell(1,superminResolution(2));
parfor i = 1:superminResolution(2)
    superminIndexArray{i} = superminArray(y(i),'y',position,...
        segmentIndex,relativeStretching,superminDistance,...
        showPlotSuperminLine);
    if verbose
        progressBar.increment(i) %#ok<PFBNS>
    end
end

if verbose
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

function superminIndexArray = superminArray(xOrY,xOrYFlag,position,...
    segmentIndex,relativeStretching,superminDistance,showPlotSuperminLine)

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

if ~isempty(intersectionPosition)
    superminIndexLine = superminimize(intersectionPosition,...
        relativeStretchingLine,superminDistance,showPlotSuperminLine);
else
    superminIndexLine = [];
end

nSupermin = length(superminIndexLine);

% Highlight superminimzed segments
if showPlotSuperminLine
    iSuperminStrainline = indexLine(1,superminIndexLine);
    iSuperminSegment = indexLine(2,superminIndexLine);
    
    superminStartIndex = arrayfun(@(idx)segmentIndex{iSuperminStrainline(idx)}(iSuperminSegment(idx),1),1:nSupermin);
    superminEndIndex = arrayfun(@(idx)segmentIndex{iSuperminStrainline(idx)}(iSuperminSegment(idx),2),1:nSupermin);
    superminPosition = arrayfun(...
        @(idx)position{iSuperminStrainline(idx)}...
        (superminStartIndex(idx):superminEndIndex(idx),:),1:nSupermin,...
        'UniformOutput',false);
    
    hAxes = showPlotSuperminLine;
    if strcmp(xOrYFlag,'x')
        hSuperminLine = plot(hAxes,xOrY*ones(1,2),get(hAxes,'ylim'));
    else
        hSuperminLine = plot(hAxes,get(hAxes,'xlim'),xOrY*ones(1,2));
    end
    set(hSuperminLine,'linestyle','--')
    
    hStrainline = cellfun(@(position)plot(hAxes,position(:,1),...
        position(:,2)),superminPosition);
    set(hStrainline,'color','r')
    pause

    if all(ishandle(hStrainline))
        delete(hStrainline)
    end
    if ishandle(hSuperminLine)
        delete(hSuperminLine)
    end
end

for i = 1:nSupermin
    superminIndexArray{indexLine(1,superminIndexLine(i))}...
        (indexLine(2,superminIndexLine(i))) = true;
end

function side = side_point_line(point,line)
%SIDE_POINT_LINE Determine on which side of a line a point lies

cellPoint = mat2cell(point,ones(size(point,1),1));

% Formula based on determinant
% Reference: http://stackoverflow.com/questions/1560492/
side = sign(cellfun(@(iPoint) det([line(2,:)-line(1,:);iPoint-line(1,:)]),cellPoint));
