function strainline_stats(strainline)

minGeodesicDeviation = min(cellfun(@min,strainline.geodesicDeviation));
maxGeodesicDeviation = max(cellfun(@max,strainline.geodesicDeviation));

sumNumel = sum(cellfun(@numel,strainline.geodesicDeviation));

% Array containing the geodesic deviation at every point.
geodesicDeviationArray = nan(1,sumNumel);

counter = 1;
for m = 1:size(strainline.geodesicDeviation,2)
    nPoints = size(strainline.geodesicDeviation{m},1);
    geodesicDeviationArray(counter:counter+nPoints-1) = ...
        strainline.geodesicDeviation{m};
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

fprintf('Geodesic deviation statistics:\n')
fprintf('\tminimum = %g\n',minGeodesicDeviation)
fprintf('\tmaximum = %g\n',maxGeodesicDeviation)
fprintf('\tmean = %g\n',meanGeodesicDeviation)
fprintf('\tmedian = %g\n',medianGeodesicDeviation)
