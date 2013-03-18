% angle_change Angular change of vector field between adjacent grid points.
%
% SYNTAX
% [thetaX,thetaY] = angle_change(vector)
% [thetaX,thetaY,thetaMax] = angle_change(vector)
%
% EXAMPLE
% bickleyJet = bickley_jet_coupled(3);
% fullDomain = bickleyJet.flow.domain;
% domain = [2.6,2.7;-.04,.01];
% bickleyJet.flow = set_flow_domain(domain,bickleyJet.flow);
% bickleyJet.flow = set_flow_resolution(64,bickleyJet.flow);
% bickleyJet.strainline = set_strainline_resolution([5 4],bickleyJet.strainline);
% showPlot.quiver = true;
% showPlot.strainline = true;
% [bickleyJet,hAxes] = strain_lcs_script(bickleyJet,showPlot);
% set(findobj(hAxes,'tag','strainline'),'color','k')
% xi1 = bickleyJet.flow.cgEigenvector(:,[1,2]);
% xi1 = shiftdim(reshape(xi1,[fliplr(bickleyJet.flow.resolution),2]),2);
% [thetaX,thetaY,thetaMax] = angle_change(xi1);
%
% print_theta_hist([thetaX(:);thetaY(:)])
%
% deltaX = diff(bickleyJet.flow.domain(1,:))/double(bickleyJet.flow.resolution(1));
% deltaY = diff(bickleyJet.flow.domain(2,:))/double(bickleyJet.flow.resolution(2));
% x = linspace(bickleyJet.flow.domain(1,1)+.5*deltaX,bickleyJet.flow.domain(1,2)-.5*deltaX,bickleyJet.flow.resolution(1)-1);
% y = linspace(bickleyJet.flow.domain(2,1)+.5*deltaY,bickleyJet.flow.domain(2,2)-.5*deltaY,bickleyJet.flow.resolution(2)-1);
% hFigure = figure;
% hAxes = axes;
% set(hAxes,'parent',hFigure)
% set(hAxes,'nextplot','add')
% set(hAxes,'DataAspectRatioMode','manual')
% set(hAxes,'DataAspectRatio',[1,1,1])
% set(hAxes,'xlim',domain(1,:))
% set(hAxes,'ylim',domain(2,:))
% hImagesc = imagesc(x,y,radtodeg(thetaMax));
% set(hImagesc,'parent',hAxes)
% hColorbar = colorbar;
% set(hColorbar,'parent',hFigure)
%
% SEE ALSO
% print_theta_hist

function [thetaX,thetaY,varargout] = angle_change(vector)

% Compare x to x + DeltaX
a = vector(:,:,1:end-1);
b = vector(:,:,2:end);

normA = norm_array(a);
normB = norm_array(b);

% Formula used is: dot(a,b) = norm(a)*norm(b)*cos(theta)
thetaX = acos(squeeze(dot(a,b,1))./normA./normB);

% Compare y to y + DeltaY
a = vector(:,1:end-1,:);
b = vector(:,2:end,:);

normA = norm_array(a);
normB = norm_array(b);

thetaY = acos(squeeze(dot(a,b,1))./normA./normB);

if nargout == 3
    % Calculate the maximum angular change at a point.
    thetaMax = nan(size(thetaY,1),size(thetaX,2));
    for m = 1:size(thetaY,1)
        for n = 1:size(thetaX,2)
            thetaMax(m,n) = max([thetaX(m,n),thetaX(m+1,n),thetaY(m,n),thetaY(m,n+1)]);
        end
    end
    varargout{1} = thetaMax;
end

% Norm of 2-by-m-by-n array.
function normArray = norm_array(array)

normArray = squeeze(sqrt(array(1,:,:).^2 + array(2,:,:).^2));
