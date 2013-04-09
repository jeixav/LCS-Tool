function shearline = geodesic_deviation_shearline(flow,shearline)

% Reshape to meshgrid array
cgEigenvalue = shiftdim(reshape(flow.cgEigenvalue,[fliplr(flow.resolution),2]),2);
cgEigenvector = shiftdim(reshape(flow.cgEigenvector,[fliplr(flow.resolution),2,2]),2);

p = inputParser;
p.KeepUnmatched = true;
addParamValue(p,'imposeIncompressibility',false,@islogical)
parse(p,flow)

[shearline.geodesicDeviationPosMeshgrid,shearline.geodesicDeviationNegMeshgrid] = shear_geodesic_deviation(cgEigenvector,cgEigenvalue,flow.domain,flow.resolution,p.Results.imposeIncompressibility);

x = linspace(flow.domain(1,1),flow.domain(1,2),flow.resolution(1));
y = linspace(flow.domain(2,1),flow.domain(2,2),flow.resolution(2));

dPosInterpolant = griddedInterpolant({y,x},shearline.geodesicDeviationPosMeshgrid);
shearline.geodesicDeviationPos = cellfun(@(position)dPosInterpolant(position(:,2),position(:,1)),shearline.positionPos,'UniformOutput',false);
shearline.averageGeodesicDeviationPos = cellfun(@average_geodesic_deviation_individual,shearline.positionPos,shearline.geodesicDeviationPos);

dNegInterpolant = griddedInterpolant({y,x},shearline.geodesicDeviationNegMeshgrid);
shearline.geodesicDeviationNeg = cellfun(@(position)dNegInterpolant(position(:,2),position(:,1)),shearline.positionNeg,'UniformOutput',false);
shearline.averageGeodesicDeviationNeg = cellfun(@average_geodesic_deviation_individual,shearline.positionNeg,shearline.geodesicDeviationNeg);

function averageGeodesicDeviation = average_geodesic_deviation_individual(shearlinePosition,geodesicDeviation)

if size(shearlinePosition,1) == 1
    averageGeodesicDeviation = geodesicDeviation;
else
    dPosition = diff(shearlinePosition);
    arcLength = [0; cumsum(hypot(dPosition(:,1),dPosition(:,2)))];
    averageGeodesicDeviation = trapz(arcLength,geodesicDeviation)/...
        arcLength(end);
end
