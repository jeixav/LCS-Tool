function [minGeodesicDeviation,maxGeodesicDeviation, ...
    meanGeodesicDeviation, medianGeodesicDeviation] = ...
    geodesic_deviation_stats(geodesicDeviation,verbose)

if nargin < 2
    verbose = false;
end

minGeodesicDeviation = min(cellfun(@min,geodesicDeviation));
maxGeodesicDeviation = max(cellfun(@max,geodesicDeviation));

sumNumel = sum(cellfun(@numel,geodesicDeviation));

% Array containing the geodesic deviation at every point.
geodesicDeviationArray = nan(1,sumNumel);

counter = 1;
for m = 1:size(geodesicDeviation,2)
    nPoints = size(geodesicDeviation{m},1);
    geodesicDeviationArray(counter:counter+nPoints-1) = ...
        geodesicDeviation{m};
    counter = counter+nPoints;
end

if sum(isnan(geodesicDeviationArray))
    nNanPoints = sum(isnan(geodesicDeviationArray));
    warning([mfilename,':GeodesicDeviationNan'],[num2str(nNanPoints),...
        ' NaN points. Total number of points: ',...
        num2str(numel(geodesicDeviationArray))])
end

meanGeodesicDeviation = ...
    mean(geodesicDeviationArray(~isnan(geodesicDeviationArray)));
medianGeodesicDeviation = ...
    median(geodesicDeviationArray(~isnan(geodesicDeviationArray)));

if verbose
    fprintf('Geodesic deviation statistics:\n')
    fprintf('\tminimum = %g\n',minGeodesicDeviation)
    fprintf('\tmaximum = %g\n',maxGeodesicDeviation)
    fprintf('\tmean = %g\n',meanGeodesicDeviation)
    fprintf('\tmedian = %g\n',medianGeodesicDeviation)
end
