function geodesicDeviation = geodesic_deviation_strainline(...
    strainlinePosition,gridPosition,maxEigenvalue,cgEigenvectors,...
    resolution)

interpolationMethod = 'linear';

l2 = reshape(maxEigenvalue,fliplr(resolution));
xi = reshape(gridPosition(:,1),fliplr(resolution));
yi = reshape(gridPosition(:,2),fliplr(resolution));
meshgridPosition = cat(3,xi,yi);
[dl2dx dl2dy] = gradient(l2,xi',yi);

v1(:,:,1) = reshape(cgEigenvectors(:,1),fliplr(resolution));
v1(:,:,2) = reshape(cgEigenvectors(:,2),fliplr(resolution));
v2(:,:,1) = reshape(cgEigenvectors(:,3),fliplr(resolution));
v2(:,:,2) = reshape(cgEigenvectors(:,4),fliplr(resolution));
[dv1xdx dv1xdy] = gradient(v1(:,:,1),xi',yi);
[dv1ydx dv1ydy] = gradient(v1(:,:,2),xi',yi);

nStrainline = size(strainlinePosition,2);
progressBar = ParforProgressStarter2(mfilename,nStrainline);

geodesicDeviation = cell(1,nStrainline);
parfor i = 1:nStrainline
    geodesicDeviation{i} = geodesic_deviation_interp(...
        strainlinePosition{i},meshgridPosition,l2,dl2dx,dl2dy,v1,v2,...
        dv1xdx,dv1xdy,dv1ydx,dv1ydy,interpolationMethod);
    progressBar.increment(i) %#ok<PFBNS>
end

try
    delete(progressBar);
catch me %#ok<NASGU>
end

end

function geodesicDeviation = geodesic_deviation_interp(interpPosition,...
    meshgridPosition,l2,dl2dx,dl2dy,v1,v2,dv1xdx,dv1xdy,dv1ydx,dv1ydy,...
    interpolationMethod)

l2_i = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),l2,...
    interpPosition(:,1),interpPosition(:,2),interpolationMethod);

dl2dx_i = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    dl2dx,interpPosition(:,1),interpPosition(:,2),interpolationMethod);
dl2dy_i = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    dl2dy,interpPosition(:,1),interpPosition(:,2),interpolationMethod);

v1x = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    v1(:,:,1),interpPosition(:,1),interpPosition(:,2),interpolationMethod);
v1y = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    v1(:,:,2),interpPosition(:,1),interpPosition(:,2),interpolationMethod);

v2x = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    v2(:,:,1),interpPosition(:,1),interpPosition(:,2),interpolationMethod);
v2y = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    v2(:,:,2),interpPosition(:,1),interpPosition(:,2),interpolationMethod);

dv1xdx_i = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    dv1xdx,interpPosition(:,1),interpPosition(:,2),interpolationMethod);
dv1xdy_i = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    dv1xdy,interpPosition(:,1),interpPosition(:,2),interpolationMethod);
dv1ydx_i = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    dv1ydx,interpPosition(:,1),interpPosition(:,2),interpolationMethod);
dv1ydy_i = interp2(meshgridPosition(:,:,1),meshgridPosition(:,:,2),...
    dv1ydy,interpPosition(:,1),interpPosition(:,2),interpolationMethod);

kappa_1 = nan(size(interpPosition,1),1);
for m = 1:size(interpPosition,1)
    grad_v1 = [dv1xdx_i(m) dv1xdy_i(m); dv1ydx_i(m) dv1ydy_i(m)];
    kappa_1(m) = dot(grad_v1*[v1x(m); v1y(m)],[v2x(m); v2y(m)]);
end

% Equation (33) from Haller and Beron-Vera, "Geodesic Theory of Transport
% Barriers in Two-Dimensional Flows", Physica D.
% FIXME Should use the general (not the incompressible) formula
geodesicDeviation = ...
    abs((dl2dx_i.*v2x + dl2dy_i.*v2y) + 2*l2_i.*kappa_1)./(2*l2_i.^3);

end
