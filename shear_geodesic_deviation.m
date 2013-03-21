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
% cgEigenvalue = shiftdim(reshape(cgEigenvalue,[fliplr(resolution),2]),2);
% cgEigenvector = shiftdim(reshape(cgEigenvector,[fliplr(resolution),2,2]),2);
% 
% [dPos,dNeg] = shear_geodesic_deviation(cgEigenvector,cgEigenvalue,domain,resolution,incompressible);
%
% Plot
% hFigure = figure;
% 
% hAxes = axes;
% set(hAxes,'parent',hFigure)
% set(hAxes,'nextplot','add')
% set(hAxes,'DataAspectRatioMode','manual')
% set(hAxes,'DataAspectRatio',[1,1,1])
% set(hAxes,'xlim',domain(1,:))
% set(hAxes,'ylim',domain(2,:))
%
% hImagesc = imagesc(domain(1,:),domain(2,:),log10(geodesicDeviation));
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
% 
% % Eta quiver plot
% n = eta(cgEigenvector,cgEigenvalue);
% x = linspace(domain(1,1),domain(1,2),resolution(1));
% y = linspace(domain(2,1),domain(2,2),resolution(2));
% hQuiver = quiver(x,y,squeeze(n(1,:,:)),squeeze(n(2,:,:)));
% set(hQuiver,'parent',hAxes)
% set(hQuiver,'AutoScaleFactor',.5)
% set(hQuiver,'color','w')

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
kappa(1,:,:) = kappa_array(squeeze(cgEigenvector(:,1,:,:)),squeeze(cgEigenvector(:,2,:,:)),gridSpacing);
kappa(2,:,:) = kappa_array(squeeze(cgEigenvector(:,2,:,:)),squeeze(cgEigenvector(:,1,:,:)),gridSpacing);

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

% FIXME Too much copy-pasting used to write this function! Will be hard to
% debug.
function k = kappa_array(cgEigenvector1,cgEigenvector2,spacing)

% Validate answers with MATLAB's gradient function
% [gradCgEigenvector2(1,1,:,:),gradCgEigenvector2(1,2,:,:)] = gradient(squeeze(cgEigenvector2(1,:,:)),spacing(1),spacing(2));
% [gradCgEigenvector2(2,1,:,:),gradCgEigenvector2(2,2,:,:)] = gradient(squeeze(cgEigenvector2(2,:,:)),spacing(1),spacing(2));

% Gradient calculation
resolution(1) = size(cgEigenvector1,3);
resolution(2) = size(cgEigenvector1,2);

gradCgEigenvector1 = nan(2,2,resolution(2),resolution(1));

% Centered finite-difference inside domain
for m = 2:resolution(1)-1
    for n = 2:resolution(2)-1
        % Reference vector to detect orientation discontinuities
        a = squeeze(cgEigenvector1(:,n,m));

        % x-direction gradient
        b = squeeze(cgEigenvector1(:,n,m+1));
        c = squeeze(cgEigenvector1(:,n,m-1));
        % FIXME For high accuracy, should not use uniform spacing
        gradCgEigenvector1(:,1,n,m) = centered_grad([c,a,b],spacing(1));
        % y-direction gradient
        b = squeeze(cgEigenvector1(:,n+1,m));
        c = squeeze(cgEigenvector1(:,n-1,m));
        gradCgEigenvector1(:,2,n,m) = centered_grad([c,a,b],spacing(2));
    end
end

% One-sided finite-difference on boundaries
for m = [1,resolution(1)]
    for n = 2:resolution(2)-1
        a = squeeze(cgEigenvector1(:,n,m));

        if m == 1
            b = squeeze(cgEigenvector1(:,n,m+1));
        elseif m == resolution(1)
            b = squeeze(cgEigenvector1(:,n,m-1));
        end
        
        theta = acos(dot(a,b)/norm(a)/norm(b));
        if theta > pi/2
            b = -b;
        end
        
        if m == 1
            gradCgEigenvector1(1,1,n,m) = (b(1)-a(1))/spacing(1);
            gradCgEigenvector1(2,1,n,m) = (b(2)-a(2))/spacing(1);
        elseif m == resolution(1)
            gradCgEigenvector1(1,1,n,m) = (a(1)-b(1))/spacing(1);
            gradCgEigenvector1(2,1,n,m) = (a(2)-b(2))/spacing(1);
        end
            
        % y-direction gradient same as non-boundary grid points
        b = squeeze(cgEigenvector1(:,n+1,m));
        c = squeeze(cgEigenvector1(:,n-1,m));
        gradCgEigenvector1(:,2,n,m) = centered_grad([c,a,b],spacing(2));
    end
end

for n = [1,resolution(2)]
    for m = 2:resolution(1)-1
        a = squeeze(cgEigenvector1(:,n,m));
        
        % x-direction gradient same as non-boundary grid points
        b = squeeze(cgEigenvector1(:,n,m+1));
        c = squeeze(cgEigenvector1(:,n,m-1));
        gradCgEigenvector1(:,1,n,m) = centered_grad([c,a,b],spacing(1));
        
        if n == 1
            b = squeeze(cgEigenvector1(:,n+1,m));
        elseif n == resolution(2)
            b = squeeze(cgEigenvector1(:,n-1,m));
        end
        
        theta = acos(dot(a,b)/norm(a)/norm(b));
        if theta > pi/2
            b = -b;
        end
        
        if n == 1
            gradCgEigenvector1(1,2,n,m) = (b(1)-a(1))/spacing(2);
            gradCgEigenvector1(2,2,n,m) = (b(2)-a(2))/spacing(2);
        elseif n == resolution(2)
            gradCgEigenvector1(1,2,n,m) = (a(1)-b(1))/spacing(2);
            gradCgEigenvector1(2,2,n,m) = (a(2)-b(2))/spacing(2);
        end

    end
end

% One-sided differences at corners
% Corner(1,1)
m = 1;
n = 1;

a = squeeze(cgEigenvector1(:,n,m));

b = squeeze(cgEigenvector1(:,n,m+1));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (b(1)-a(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (b(2)-a(2))/spacing(1);

b = squeeze(cgEigenvector1(:,n+1,m));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (b(1)-a(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (b(2)-a(2))/spacing(2);

% Corner(1,res(2))
n = resolution(2);

a = squeeze(cgEigenvector1(:,n,m));

b = squeeze(cgEigenvector1(:,n,m+1));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (b(1)-a(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (b(2)-a(2))/spacing(1);

b = squeeze(cgEigenvector1(:,n-1,m));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (a(1)-b(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (a(2)-b(2))/spacing(2);

% Corner(res(1),res(2))
m = resolution(1);

a = squeeze(cgEigenvector1(:,n,m));

b = squeeze(cgEigenvector1(:,n,m-1));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (a(1)-b(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (a(2)-b(2))/spacing(1);

b = squeeze(cgEigenvector1(:,n-1,m));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (a(1)-b(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (a(2)-b(2))/spacing(2);

% Corner(res(1),1)
n = 1;

a = squeeze(cgEigenvector1(:,n,m));

b = squeeze(cgEigenvector1(:,n,m-1));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (a(1)-b(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (a(2)-b(2))/spacing(1);

b = squeeze(cgEigenvector1(:,n+1,m));

theta = acos(dot(a,b)/norm(a)/norm(b));
if theta > pi/2
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (b(1)-a(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (b(2)-a(2))/spacing(2);

% Equation 30 from doi:10.1016/j.physd.2012.06.012
k = (squeeze(gradCgEigenvector1(1,1,:,:)).*squeeze(cgEigenvector1(1,:,:)) + squeeze(gradCgEigenvector1(1,2,:,:)).*squeeze(cgEigenvector1(2,:,:))).*squeeze(cgEigenvector2(1,:,:)) + (squeeze(gradCgEigenvector1(2,1,:,:)).*squeeze(cgEigenvector1(1,:,:)) + squeeze(gradCgEigenvector1(2,2,:,:)).*squeeze(cgEigenvector1(2,:,:))).*squeeze(cgEigenvector2(2,:,:));

function grad = centered_grad(v,spacing)

theta = acos(dot(v(:,2),v(:,3))/norm(v(:,2))/norm(v(:,3)));
if theta > pi/2
    v(:,3) = -v(:,3);
end

theta = acos(dot(v(:,2),v(:,1))/norm(v(:,2))/norm(v(:,1)));
if theta > pi/2
    v(:,1) = -v(:,1);
end
grad(1) = .5*(v(1,3)-v(1,1))/spacing;
grad(2) = .5*(v(2,3)-v(2,1))/spacing;
