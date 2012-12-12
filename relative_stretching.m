%RELATIVE_STRETCHING Relative stretching of a strainline segment.
%   The relative stretching "metric" is:
%   currentLength/initialLength

function relativeStretching = relative_stretching(strainlinePosition,...
    segmentIndex,position,maxEigenvalue,resolution,verbose)

relativeStretching = cell(size(segmentIndex));

nStrainlines = size(strainlinePosition,2);

if verbose
    if ~exist('ParforProgressStarter2','file')
        addpath('ParforProgress2')
    end
    progressBar = ParforProgressStarter2(mfilename,nStrainlines);
else
    progressBar = [];
end

parfor iStrainline = 1:nStrainlines

    if ~isempty(segmentIndex{iStrainline})
        relativeStretching{iStrainline} = ...
            nan(size(segmentIndex{iStrainline},1),1);
        for iSegment = 1:size(segmentIndex{iStrainline},1)
            initialPosition = [...
                strainlinePosition{iStrainline}(segmentIndex{iStrainline}(iSegment,1):segmentIndex{iStrainline}(iSegment,2),1),...
                strainlinePosition{iStrainline}(segmentIndex{iStrainline}(iSegment,1):segmentIndex{iStrainline}(iSegment,2),2)];

            xD = diff(initialPosition(:,1));
            yD = diff(initialPosition(:,2));
            arcLength = [0; cumsum(sqrt(xD.^2 + yD.^2))];
            lambda_1 = interp2(...
                reshape(position(:,1),fliplr(resolution)),...
                reshape(position(:,2),fliplr(resolution)),...
                reshape(maxEigenvalue,fliplr(resolution)),...
                initialPosition(:,1),initialPosition(:,2)); %#ok<PFBNS>
            
            % Eliminate repeated points at beginning and end of strainline
            % FIXME These repeated points should not be generated
            % in the first place
            startIndex = find(arcLength > 0,1) - 1;
            arcLength = arcLength(startIndex:end);
            lambda_1 = lambda_1(startIndex:end);
            endIndex = find(~diff(arcLength),1);
            if endIndex
                arcLength = arcLength(1:endIndex);
                lambda_1 = lambda_1(1:endIndex);
            end
            currentLength = trapz(arcLength,sqrt(lambda_1));
            relativeStretching{iStrainline}(iSegment)...
                = currentLength/arcLength(end);
        end
    end
    if verbose
        progressBar.increment(iStrainline) %#ok<PFBNS>
    end
end

if verbose
    try
        delete(progressBar);
    catch me %#ok<NASGU>
    end
end
