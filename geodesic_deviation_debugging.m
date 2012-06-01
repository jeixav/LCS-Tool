function geodesic_deviation_debugging(flow,shearline)

idx = 1;

dPosition = diff(shearline.positionPos{idx});
arcLength = [0; cumsum(hypot(dPosition(:,1),dPosition(:,2)))];

figure
a1 = axes('nextplot','add','box','on');
plot(a1,arcLength,shearline.geodesicDeviationPos{idx},'-o')

% idx2 = [40:65 165:190];
% plot(a1,arcLength(idx2),doubleGyre.shearline.geodesicDeviationPos{idx}...
%     (idx2),'ro','markerFaceColor','r')

% figure(1)
% p1 = plot(doubleGyre.shearline.positionPos{idx}(idx2,1),...
%     doubleGyre.shearline.positionPos{idx}(idx2,2),'ro','markerFaceColor','r');

% downsampling = 10;
% dsResolution = flow.resolution/downsampling;
% 
% gridPosition.x = linspace(flow.domain(1,1),flow.domain(1,2),... 
%     dsResolution(1));
% gridPosition.y = linspace(flow.domain(2,1),flow.domain(2,2),...
%     dsResolution(2));
% 
% a1 = get(figure(1),'children');
% 
% xi1(:,:,1) = reshape(flow.cgEigenvector(:,1),fliplr(flow.resolution));
% xi1(:,:,2) = reshape(flow.cgEigenvector(:,2),fliplr(flow.resolution));
% 
% scale = .5;
% quiver(a1,gridPosition.x,gridPosition.y,...
%     xi1(1:downsampling:end,1:downsampling:end,1),...
%     xi1(1:downsampling:end,1:downsampling:end,2),...
%     scale,'tag','xi1Quiver')
% 
% axis(a1,'tight','equal')
% set(a1,'xLim',[0 1],'yLim',[.3 .85])
