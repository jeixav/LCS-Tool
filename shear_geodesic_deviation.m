% shear_geodesic_deviation Shear geodesic deviation
%
% SYNTAX
% [dPos,dNeg] = shear_geodesic_deviation(cgEigenvector,cgEigenvalue,domain,resolution,incompressible)
%
% EXAMPLE
% bickleyJet = bickley_jet(3);
% [cgEigenvalue,cgEigenvector] = eig_cgStrain(bickleyJet.flow);
%
% % Reshape to m-by-n array
% cgEigenvalue = shiftdim(reshape(cgEigenvalue,[fliplr(bickleyJet.flow.resolution),2]),2);
% cgEigenvector = shiftdim(reshape(cgEigenvector,[fliplr(bickleyJet.flow.resolution),2,2]),2);
% 
% [dPos,dNeg] = shear_geodesic_deviation(cgEigenvector,cgEigenvalue,bickleyJet.flow.domain,bickleyJet.flow.resolution,bickleyJet.flow.incompressible);
%
% % Plot
% hFigure = figure;
% 
% hAxes = axes;
% set(hAxes,'parent',hFigure)
% set(hAxes,'nextplot','add')
% set(hAxes,'DataAspectRatioMode','manual')
% set(hAxes,'DataAspectRatio',[1,1,1])
% set(hAxes,'xlim',bickleyJet.flow.domain(1,:))
% set(hAxes,'ylim',bickleyJet.flow.domain(2,:))
%
% hImagesc = imagesc(bickleyJet.flow.domain(1,:),bickleyJet.flow.domain(2,:),log10(dPos));
% set(hImagesc,'parent',hAxes);
% uistack(hImagesc,'bottom')
% 
% hXlabel = xlabel('x');
% set(hXlabel,'parent',hAxes)
% 
% hYlabel = ylabel('y');
% set(hYlabel,'parent',hAxes)
% 
% cbar_axes = colorbar;
% set(cbar_axes,'parent',get(hAxes,'parent'))
% set(get(cbar_axes,'xlabel'),'string','log(d_g)')

function [dPos,dNeg] = shear_geodesic_deviation(cgEigenvector,cgEigenvalue,domain,resolution,incompressible)

% Equations from doi:10.1016/j.physd.2012.06.012 page 1685
alpha = squeeze(sqrt(sqrt(cgEigenvalue(2,:,:))./(sqrt(cgEigenvalue(1,:,:)) + sqrt(cgEigenvalue(2,:,:)))));
beta = squeeze(sqrt(sqrt(cgEigenvalue(1,:,:))./(sqrt(cgEigenvalue(1,:,:)) + sqrt(cgEigenvalue(2,:,:)))));

xGrid = linspace(domain(1,1),domain(1,2),resolution(1));
yGrid = linspace(domain(2,1),domain(2,2),resolution(2));

[gradAlpha(1,:,:),gradAlpha(2,:,:)] = gradient(alpha,xGrid,yGrid);

% Equation 18 from doi:10.1016/j.physd.2012.06.012
etaPos(1,:,:) = alpha.*squeeze(cgEigenvector(1,1,:,:)) + beta.*squeeze(cgEigenvector(2,1,:,:));
etaPos(2,:,:) = alpha.*squeeze(cgEigenvector(1,2,:,:)) + beta.*squeeze(cgEigenvector(2,2,:,:));
etaNeg(1,:,:) = alpha.*squeeze(cgEigenvector(1,1,:,:)) - beta.*squeeze(cgEigenvector(2,1,:,:));
etaNeg(2,:,:) = alpha.*squeeze(cgEigenvector(1,2,:,:)) - beta.*squeeze(cgEigenvector(2,2,:,:));

dotGradAlphaEtaPos = squeeze(gradAlpha(1,:,:).*etaPos(1,:,:) + gradAlpha(2,:,:).*etaPos(2,:,:));
dotGradAlphaEtaNeg = squeeze(gradAlpha(1,:,:).*etaNeg(1,:,:) + gradAlpha(2,:,:).*etaNeg(2,:,:));

% FIXME Should not use uniform grid spacing
gridSpacing = [mean(diff(xGrid)),mean(diff(yGrid))];
kappa(1,:,:) = euclidean_curvature(squeeze(cgEigenvector(:,1,:,:)),squeeze(cgEigenvector(:,2,:,:)),gridSpacing);
kappa(2,:,:) = euclidean_curvature(squeeze(cgEigenvector(:,2,:,:)),squeeze(cgEigenvector(:,1,:,:)),gridSpacing);

if incompressible
    [gradCgEigenvalue2(1,:,:),gradCgEigenvalue2(2,:,:)] = gradient(squeeze(cgEigenvalue(2,:,:)),xGrid,yGrid);
    dotGradCgEigenvalue2CgEigenvector1 = squeeze(gradCgEigenvalue2(1,:,:)).*squeeze(cgEigenvector(1,1,:,:)) + squeeze(gradCgEigenvalue2(2,:,:)).*squeeze(cgEigenvector(1,2,:,:));
    dotGradCgEigenvalue2CgEigenvector2 = squeeze(gradCgEigenvalue2(1,:,:)).*squeeze(cgEigenvector(2,1,:,:)) + squeeze(gradCgEigenvalue2(2,:,:)).*squeeze(cgEigenvector(2,2,:,:));
    
    % Equation 34 from doi:10.1016/j.physd.2012.06.012
    dPos = squeeze((sqrt(1 + cgEigenvalue(2,:,:)) - sqrt(cgEigenvalue(2,:,:)))./sqrt(1 + cgEigenvalue(2,:,:))) + abs(.5*dotGradCgEigenvalue2CgEigenvector1./squeeze(cgEigenvalue(2,:,:).*sqrt(1 + cgEigenvalue(2,:,:))) - .5*dotGradCgEigenvalue2CgEigenvector2.*squeeze(((1 + cgEigenvalue(2,:,:)).^1.5 - cgEigenvalue(2,:,:).^2.5)./(cgEigenvalue(2,:,:).^3.*(1 + cgEigenvalue(2,:,:)).^1.5)) - squeeze(kappa(1,:,:)).*squeeze((cgEigenvalue(2,:,:).^2.5 + (1 - cgEigenvalue(2,:,:).^2).*sqrt(1 + cgEigenvalue(2,:,:)))./(cgEigenvalue(2,:,:).^2.*sqrt(1 + cgEigenvalue(2,:,:)))) + squeeze(kappa(2,:,:))./squeeze(sqrt(1 + cgEigenvalue(2,:,:))));
    dNeg = squeeze((sqrt(1 + cgEigenvalue(2,:,:)) - sqrt(cgEigenvalue(2,:,:)))./sqrt(1 + cgEigenvalue(2,:,:))) + abs(.5*dotGradCgEigenvalue2CgEigenvector1./squeeze(cgEigenvalue(2,:,:).*sqrt(1 + cgEigenvalue(2,:,:))) + .5*dotGradCgEigenvalue2CgEigenvector2.*squeeze(((1 + cgEigenvalue(2,:,:)).^1.5 - cgEigenvalue(2,:,:).^2.5)./(cgEigenvalue(2,:,:).^3.*(1 + cgEigenvalue(2,:,:)).^1.5)) + squeeze(kappa(1,:,:)).*squeeze((cgEigenvalue(2,:,:).^2.5 + (1 - cgEigenvalue(2,:,:).^2).*sqrt(1 + cgEigenvalue(2,:,:)))./(cgEigenvalue(2,:,:).^2.*sqrt(1 + cgEigenvalue(2,:,:)))) + squeeze(kappa(2,:,:))./squeeze(sqrt(1 + cgEigenvalue(2,:,:))));
else
    [gradCgEigenvalue1(1,:,:),gradCgEigenvalue1(2,:,:)] = gradient(squeeze(cgEigenvalue(1,:,:)),xGrid,yGrid);
    dotGradCgEigenvalue1CgEigenvector2 = squeeze(gradCgEigenvalue1(1,:,:)).*squeeze(cgEigenvector(1,2,:,:)) + squeeze(gradCgEigenvalue1(2,:,:)).*squeeze(cgEigenvector(2,2,:,:));
    
    % Equation 33 from doi:10.1016/j.physd.2012.06.012
    dPos = abs(1-alpha) + abs(-dotGradAlphaEtaPos./beta + alpha.*squeeze(kappa(1,:,:)) - beta.*squeeze(kappa(2,:,:)) + (squeeze(cgEigenvalue(1,:,:))./squeeze(cgEigenvalue(2,:,:)) - 1).*squeeze(kappa(1,:,:)) - .5*dotGradCgEigenvalue1CgEigenvector2./squeeze(cgEigenvalue(2,:,:)));
    dNeg = abs(1-alpha) + abs( dotGradAlphaEtaNeg./beta + alpha.*squeeze(kappa(1,:,:)) + beta.*squeeze(kappa(2,:,:)) + (squeeze(cgEigenvalue(1,:,:))./squeeze(cgEigenvalue(2,:,:)) - 1).*squeeze(kappa(1,:,:)) - .5*dotGradCgEigenvalue1CgEigenvector2./squeeze(cgEigenvalue(2,:,:)));
end

% euclidean_curvature Euclidean curvature of a trajectory of the vector2
% vector field
%
% SYNTAX
% kappa = euclidean_curvature(vector1,vector2,spacing)
%
% DESCRIPTION
% Equation (30) from doi:10.1016/j.physd.2012.06.012. Orientation
% discontinuities in the eigenvector field are taken into account.

function kappa = euclidean_curvature(vector1,vector2,spacing)

rotation_check(vector1,vector2)

resolution(1) = size(vector1,3);
resolution(2) = size(vector1,2);

gradCgEigenvector1 = nan(2,2,resolution(2),resolution(1));

% x-direction gradient
ref = vector1(:,:,2:end-1);

refP = vector1(:,:,3:end);
s = sign(squeeze(dot(ref,refP)));
s = shiftdim(repmat(s,[1,1,2]),2);
refP = s.*refP;

refM = vector1(:,:,1:end-2);
s = sign(squeeze(dot(ref,refM)));
s = shiftdim(repmat(s,[1,1,2]),2);
refM = s.*refM;

gradCgEigenvector1(:,1,:,2:end-1) = .5*(refP - refM)/spacing(1);

% y-direction gradient
ref = vector1(:,2:end-1,:);

refP = vector1(:,3:end,:);
s = sign(squeeze(dot(ref,refP)));
s = shiftdim(repmat(s,[1,1,2]),2);
refP = s.*refP;

refM = vector1(:,1:end-2,:);
s = sign(squeeze(dot(ref,refM)));
s = shiftdim(repmat(s,[1,1,2]),2);
refM = s.*refM;

gradCgEigenvector1(:,2,2:end-1,:) = .5*(refP - refM)/spacing(2);

% One-sided finite-difference on boundaries
for m = [1,resolution(1)]
    ref = vector1(:,:,m);

    if m == 1
        refP = vector1(:,:,m+1);
        s = sign(squeeze(dot(ref,refP)));
    else
        refM = vector1(:,:,m-1);
        s = sign(squeeze(dot(ref,refM)));
    end
    s = repmat(s,[2,1]);

    if m == 1
        refP = s.*refP;
        gradCgEigenvector1(:,1,:,m) = (refP-ref)/spacing(1);
    else
        refM = s.*refM;
        gradCgEigenvector1(:,1,:,m) = (ref-refM)/spacing(1);
    end
end

for n = [1,resolution(2)]
    ref = squeeze(vector1(:,n,:));
    
    if n == 1
        refP = squeeze(vector1(:,n+1,:));
        s = sign(squeeze(dot(ref,refP)));
    else
        refM = squeeze(vector1(:,n-1,:));
        s = sign(squeeze(dot(ref,refM)));
    end
    s = repmat(s,[2,1]);

    if n == 1
        refP = s.*refP;
        gradCgEigenvector1(:,2,n,:) = (refP-ref)/spacing(2);
    else
        refM = s.*refM;
        gradCgEigenvector1(:,2,n,:) = (ref-refM)/spacing(2);
    end
end

% One-sided differences at corners
% Corner(1,1)
m = 1;
n = 1;

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m+1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (b(1)-a(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (b(2)-a(2))/spacing(1);

b = squeeze(vector1(:,n+1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (b(1)-a(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (b(2)-a(2))/spacing(2);

% Corner(1,res(2))
n = resolution(2);

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m+1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (b(1)-a(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (b(2)-a(2))/spacing(1);

b = squeeze(vector1(:,n-1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (a(1)-b(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (a(2)-b(2))/spacing(2);

% Corner(res(1),res(2))
m = resolution(1);

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m-1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (a(1)-b(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (a(2)-b(2))/spacing(1);

b = squeeze(vector1(:,n-1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (a(1)-b(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (a(2)-b(2))/spacing(2);

% Corner(res(1),1)
n = 1;

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m-1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (a(1)-b(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (a(2)-b(2))/spacing(1);

b = squeeze(vector1(:,n+1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (b(1)-a(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (b(2)-a(2))/spacing(2);

% Equation 30 from doi:10.1016/j.physd.2012.06.012
kappa = (squeeze(gradCgEigenvector1(1,1,:,:)).*squeeze(vector1(1,:,:)) + squeeze(gradCgEigenvector1(1,2,:,:)).*squeeze(vector1(2,:,:))).*squeeze(vector2(1,:,:)) + (squeeze(gradCgEigenvector1(2,1,:,:)).*squeeze(vector1(1,:,:)) + squeeze(gradCgEigenvector1(2,2,:,:)).*squeeze(vector1(2,:,:))).*squeeze(vector2(2,:,:));

% Ensure when vector1 has an orientation discontinuity, so does vector2.
% Check this by verifying vector2 is consistently rotated in the same
% direction (by ±90°) compared to vector1.
%
% The custom method used to calculate eigenvalues should guarantee this
% always. MATLAB's eig function seems to satisfy this condition as well,
% but it is not documented, so it is best to verify.
function rotation_check(vector1,vector2)

% Set rotation based on that of first grid point vectors
rotationMatrix = [0,1;-1,0];
if any(rotationMatrix*vector2(:,1) ~= vector1(:,1))
    rotationMatrix = [0,-1;1,0];
    if any(rotationMatrix*vector2(:,1) ~= vector1(:,1))
        error(['Unable to determine 90° rotation direction from vector1(:,1) = ',num2str(transpose(vector1(:,1))),' and vector1(:,2) = ',num2str(transpose(vector2(:,1))),'.'])
    end
end

nanIdx = isnan(vector1(1,:));

if any(any(rotationMatrix*vector2(:,~nanIdx) ~= vector1(:,~nanIdx)))
    error('vector1 not rotated by 90° with respect to vector2')
end
