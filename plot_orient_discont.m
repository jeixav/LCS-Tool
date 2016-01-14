% plot_orient_discont Plot eigenvector field orientation discontinuities
%
% SYNTAX
% hAxes = plot_orient_discont(eigenvector,domain,resolution)
% hAxes = plot_orient_discont(eigenvector,domain,resolution,angleTol)
%
% EXAMPLE
%   s = load('datasets/bickley_jet/bickleyJet1.mat');
%   eigenvector = [s.bickleyJet.flow.cgEigenvector(:,1),s.bickleyJet.flow.cgEigenvector(:,2)];
%   a1 = plot_orient_discont(eigenvector,s.bickleyJet.flow.domain,s.bickleyJet.flow.resolution);
%   zoomDomain = [2.6,2.7;-.04,.01];
%   set(a1,'xlim',zoomDomain(1,:),'ylim',zoomDomain(2,:))

function hAxes = plot_orient_discont(eigenvector,domain,resolution,angleTol)

narginchk(3,4)

if nargin == 3
    angleTol = deg2rad(90);
end

f1 = figure;
hAxes = axes('parent',f1,...
    'box','on',...
    'NextPlot','add',...
    'DataAspectRatio',[1 1 1],...
    'DataAspectRatioMode','manual');

position = initialize_ic_grid(resolution,domain);
hQuiver = quiver(hAxes,position(:,1),position(:,2),eigenvector(:,1),eigenvector(:,2));
set(hQuiver,'AutoScaleFactor',.5)

positionGridX = reshape(position(:,1),fliplr(resolution));
positionGridY = reshape(position(:,2),fliplr(resolution));

eigenvector = reshape(eigenvector,[fliplr(resolution) 2]);

idx1Discont = idx1_discont(eigenvector(:,:,1),eigenvector(:,:,2),resolution,angleTol);
deltaX = .5*diff(domain(1,:))/(double(resolution(1)) - 1);
hDiscont = plot(hAxes,positionGridX(idx1Discont)+deltaX,positionGridY(idx1Discont),'ro');
set(hDiscont,'MarkerFaceColor','r')

idx2Discont = idx2_discont(eigenvector(:,:,1),eigenvector(:,:,2),resolution,angleTol);
deltaY = .5*diff(domain(2,:))/(double(resolution(2)) - 1);
hDiscont = plot(hAxes,positionGridX(idx2Discont),positionGridY(idx2Discont)+deltaY,'ro');
set(hDiscont,'MarkerFaceColor','r')

function idx1Discont = idx1_discont(v1,v2,resolution,angleTol)

idx1 = 1:resolution(1)-1;

dotProduct = v1(:,idx1).*v1(:,idx1+1) + v2(:,idx1).*v2(:,idx1+1);

idx1Discont = find(acos(dotProduct) > (pi - angleTol));

function idx2Discont = idx2_discont(v1,v2,resolution,angleTol)

idx2 = 1:resolution(2)-1;

dotProduct = v1(idx2,:).*v1(idx2+1,:) + v2(idx2,:).*v2(idx2+1,:);

idx2Discont = find(acos(dotProduct) > (pi - angleTol));
