%RELATIVE_STRETCHING Relative stretching of a strainline segment.
%   The relative stretching "metric" is:
%   currentLength/initialLength

function relativeStretching = relative_stretching(strainlinePosition,segmentIndex,gridEigenvalue,flowDomain,flowResolution,verbose)

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

x = linspace(flowDomain(1,1),flowDomain(1,2),flowResolution(1));
y = linspace(flowDomain(2,1),flowDomain(2,2),flowResolution(2));

eigenvalue = griddedInterpolant({y,x},reshape(gridEigenvalue,fliplr(flowResolution)));

for iStrainline = 1:nStrainlines
    if ~isempty(segmentIndex{iStrainline})
        relativeStretching{iStrainline} = nan(size(segmentIndex{iStrainline},1),1);
        for iSegment = 1:size(segmentIndex{iStrainline},1)
            initialPosition = [strainlinePosition{iStrainline}(segmentIndex{iStrainline}(iSegment,1):segmentIndex{iStrainline}(iSegment,2),1),strainlinePosition{iStrainline}(segmentIndex{iStrainline}(iSegment,1):segmentIndex{iStrainline}(iSegment,2),2)];

            xD = diff(initialPosition(:,1));
            yD = diff(initialPosition(:,2));
            arcLength = [0; cumsum(sqrt(xD.^2 + yD.^2))];
            eigenvalueLocal = eigenvalue(initialPosition(:,2),initialPosition(:,1));
            
            currentLength = trapz(arcLength,sqrt(eigenvalueLocal));
            relativeStretching{iStrainline}(iSegment) = currentLength/arcLength(end);
        end
    end
    if verbose
        progressBar.increment(iStrainline)
    end
end

if verbose
    try
        delete(progressBar);
    catch me %#ok<NASGU>
    end
end
