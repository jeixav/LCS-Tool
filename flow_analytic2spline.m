function pp = flow_analytic2spline(flow,resolution)
% Create spline from analytically defined flow.

resolutionX = resolution(1);
resolutionY = resolution(2);
resolutionTime = resolution(3);

xGrid = linspace(flow.domain(1,1),flow.domain(1,2),resolutionX);
yGrid = linspace(flow.domain(2,1),flow.domain(2,2),resolutionY);
[yy,xx] = meshgrid(yGrid,xGrid);
position = [xx(:) yy(:)];

timeGrid = linspace(flow.timespan(1),flow.timespan(2),resolutionTime);

F = nan(2,resolutionX,resolutionY,resolutionTime);

for iTime = 1:length(timeGrid)
    
    time = timeGrid(iTime);
   
    f = arrayfun(@(row)flow.derivative(time,position(row,:)),...
        1:resolutionX*resolutionY,'UniformOutput',false);
    f1 = cellfun(@(x)x(1),f);
    f2 = cellfun(@(x)x(2),f);

    F(1,:,:,iTime) = reshape(f1,resolutionX,resolutionY);
    F(2,:,:,iTime) = reshape(f2,resolutionX,resolutionY);

end

pp = csapi({xGrid,yGrid,timeGrid},F);