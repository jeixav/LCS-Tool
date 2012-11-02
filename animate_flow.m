% animate_flow Display flow animation
%
% DESCRIPTION
% flow = animate_flow(flow,animationTime,framerate)
% animationTime has units of seconds
% framerate has units of 1/second
%
% EXAMPLE
% addpath('flow_templates')
% doubleGyre = double_gyre;
% doubleGyre.flow = animate_flow(doubleGyre.flow);

function flow = animate_flow(flow,animationTime,framerate)

if nargin < 3
    framerate = 10;
end

if nargin < 2
    animationTime = 10;
end

initialPosition = initialize_ic_grid(flow.resolution,flow.domain);
if isfield(flow,'symDerivative') && ~isfield(flow,'derivative')
    flow.derivative = sym2fun(flow.symDerivative);
end
solution = integrate_flow(flow,initialPosition);

mainAxes = setup_figure(flow);

position = arrayfun(@(iSolution)deval(iSolution,flow.timespan(1)),...
    solution,'UniformOutput',false);
position = cell2mat(position);
p1 = plot(mainAxes,position(1,:),position(2,:),'o');

t1 = text(0,0,['t = ',num2str(flow.timespan(1))]);
xLabelPos = get(get(mainAxes,'XLabel'),'position');
t1Pos = [flow.domain(1,1) + diff(flow.domain(1,:))*.9 xLabelPos(2)];
set(t1,'Parent',mainAxes,...
    'position',t1Pos,...
    'VerticalAlignment',get(get(mainAxes,'XLabel'),'verticalalignment'),...
    'Tag','timeText')
drawnow
tic

totalFrames = framerate*animationTime+1;
timesteps = linspace(flow.timespan(1),flow.timespan(2),totalFrames);

saveAnimation = false;
if saveAnimation
    % frame = struct(1,totalFrames);
    writerObj = VideoWriter('test_animation.avi');
    writerObj.FrameRate = framerate;
    open(writerObj)
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end

for idx = 2:length(timesteps)
    position = arrayfun(@(iSolution)deval(iSolution,timesteps(idx)),...
        solution,'UniformOutput',false);
    position = cell2mat(position);

    if isfield(flow,'periodicBc') && any(flow.periodicBc)
        position = impose_periodic_bc2(position,flow.domain,...
            flow.periodicBc);
    end

    set(p1,'xData',position(1,:),'yData',position(2,:))
    set(t1,'string',['t = ',num2str(timesteps(idx))])
    drawnow
    delay = 1/framerate - toc;
    if delay < 0
        warning('animate_flow:negative_delay','Frame rate too high')
    end
    if saveAnimation
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    pause(delay)
    tic
end

if saveAnimation
    close(writerObj)
end
