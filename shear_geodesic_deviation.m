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

