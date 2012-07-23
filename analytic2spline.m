function pp = analytic2spline(flow,resolution)
%ANALYTIC2SPLINE  Create spline from analytically defined flow.
%   PP = ANALYTIC2SPLINE(FLOW,RESOLUTION) takes the analytic flow
%   definition in FLOW, evaluates it at RESOLUTION = [RESOLUTIONX
%   RESOLUTIONY RESOLUTIONTIME] and creates the spline PP.
%
%   Example:
%   addpath('flow_templates')
%   doubleGyre = double_gyre;
%   resolutionX = 20;
%   resolutionY = 10;
%   resolutionTime = 40;
% 
%   resolution = [resolutionX resolutionY resolutionTime];
%
%   pp = analytic2spline(doubleGyre.flow,resolution);
% 
%   derivative = @(t,y)fnval(pp,[y(1); y(2); t]);
%   dDerivative = @(t,y)reshape(fnval(fndir(pp,[eye(2); 0 0]),...
%     [y(1); y(2); t]),2,2);
% 
%   norm(abs(doubleGyre.flow.derivative(time,position) ...
%     - derivative(time,position)))
%   norm(abs(doubleGyre.flow.dDerivative(time,position) ...
%     - dDerivative(time,position)))

resolutionX = resolution(1);
resolutionY = resolution(2);
resolutionTime = resolution(3);

% FIXME Unlikely to get good results if resolution = 1
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