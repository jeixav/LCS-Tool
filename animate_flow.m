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

function flow = animate_flow(flow,animationTime,framerate,...
    animationFilename)

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
if isfield(flow,'symDerivative') && ~isfield(flow,'derivative')
    flow.derivative = sym2fun(flow.symDerivative);
end

if ~isfield(flow,'solution')
    flow.solution = integrate_flow(flow,initialPosition);
end

mainAxes = setup_figure(flow.domain);

position = arrayfun(@(iSolution)deval(iSolution,flow.timespan(1)),...
    flow.solution,'UniformOutput',false);
position = cell2mat(position);
p1 = plot(mainAxes,position(1,:),position(2,:));
set(p1,'LineStyle','none')
set(p1,'Marker','o')
set(p1,'MarkerFaceColor','b')
set(p1,'MarkerSize',4)

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
    position = arrayfun(@(iSolution)deval(iSolution,timesteps(idx)),...
        flow.solution,'UniformOutput',false);
    position = cell2mat(position);

    if isfield(flow,'periodicBc') && any(flow.periodicBc)
        position = impose_periodic_bc2(position,flow.domain,...
            flow.periodicBc);
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
