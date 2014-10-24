% animate_flow Display flow animation
%
% DESCRIPTION
% flow = animate_flow(flow,animationTime,framerate,animationFilename)
% animationTime has units of seconds
% framerate has units of 1/second
%
% EXAMPLE
% addpath('flow_templates')
% doubleGyre = double_gyre;
% doubleGyre.flow = animate_flow(doubleGyre.flow);

function flow = animate_flow(flow,animationTime,framerate,animationFilename)

narginchk(1,4)

% FIXME Replace with inputParser
if nargin < 4
    animationFilename = false;
end

if nargin < 3
    framerate = 10;
end

if nargin < 2
    animationTime = 10;
end

initialPosition = initialize_ic_grid(flow.resolution,flow.domain);

if ~isfield(flow,'solution')
    % The equation of variation terms are not needed to produce the flow
    % animation
    useEoV = false;
    flow.solution = integrate_flow(flow,initialPosition,useEoV);
end

mainAxes = setup_figure(flow.domain);

% Set default values for flow structure
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'coupledIntegration',true,@islogical);
parse(p,flow);
coupledIntegration = p.Results.coupledIntegration;

if coupledIntegration
    position = deval(flow.solution,flow.timespan(1));
    position = [position(1:2:end-1),position(2:2:end)];
    position = transpose(position);
else
    position = arrayfun(@(iSolution)deval(iSolution,flow.timespan(1)),flow.solution,'UniformOutput',false);
    position = cell2mat(position);
end
p1 = plot(mainAxes,position(1,:),position(2,:));
set(p1,'LineStyle','none')
set(p1,'Marker','o')
set(p1,'MarkerFaceColor','b')
set(p1,'MarkerSize',4)

t1 = text(0,0,['t = ',num2str(flow.timespan(1))]);
xLabelPos = get(get(mainAxes,'XLabel'),'position');
t1Pos = [flow.domain(1,1) + diff(flow.domain(1,:))*.9 xLabelPos(2)];
set(t1,'parent',mainAxes)
set(t1,'position',t1Pos)
set(t1,'VerticalAlignment',get(get(mainAxes,'XLabel'),'VerticalAlignment'))
set(t1,'tag','timeText')
drawnow
tic

totalFrames = framerate*animationTime+1;
timesteps = linspace(flow.timespan(1),flow.timespan(2),totalFrames);

if animationFilename
    if exist(animationFilename,'file')
        error(['File exists: ',animationFilename])
    end
    writerObj = VideoWriter(animationFilename);
    writerObj.FrameRate = framerate;
    open(writerObj)
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

warningDisplayed = false;
for idx = 2:length(timesteps)
    if coupledIntegration
        position = deval(flow.solution,timesteps(idx));
        position = [position(1:2:end-1),position(2:2:end)];
        position = transpose(position);
    else
        position = arrayfun(@(iSolution)deval(iSolution,timesteps(idx)),flow.solution,'UniformOutput',false);
        position = cell2mat(position);
    end
    
    if isfield(flow,'periodicBc') && any(flow.periodicBc)
        position = apply_periodic_bc(position,flow.periodicBc,flow.domain);
    end

    set(p1,'xData',position(1,:),'yData',position(2,:))
    set(t1,'string',['t = ',num2str(timesteps(idx))])
    drawnow
    delay = 1/framerate - toc;
    if (delay < 0) && (warningDisplayed == false)
        warning([mfilename,':negative_delay'],'Frame rate too high')
        warningDisplayed = true;
    end
    if animationFilename
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    pause(delay)
    tic
end

if animationFilename
    close(writerObj)
end
