function flow = animate_flow(flow,verbose)

if nargin < 2
    verbose.progress = true;
end

if ~isfield(flow,'solution')
    
    if verbose.progress
        cpb = ConsoleProgressBar;
        cpb.setText(mfilename)
        cpb.setTextPosition('left')
        cpb.setElapsedTimeVisible(1);
        cpb.setRemainedTimeVisible(1);
        cpb.setLength(20)
        cpb.start
    end

    initialPosition = initial_position(flow.domain,flow.resolution);
    flow.solution = integrate_flow2(flow,initialPosition);
    
    if verbose.progress
        cpb.setValue(cpb.maximum)
        cpb.stop
        fprintf('\n')
    end

end

mainAxes = setup_figure(flow);

position = arrayfun(@(iSolution)deval(iSolution,flow.timespan(1)),...
    flow.solution,'UniformOutput',false);
position = cell2mat(position);
p1 = plot(mainAxes,position(1,:),position(2,:),'o');

t1 = text(0,0,['t = ',num2str(flow.timespan(1))]);
set(t1,'Parent',mainAxes,...
    'units','normalized',...
    'position',[.9 -.2 0],...
    'Tag','timeText')
drawnow
tic

framesPerSecond = 10;
timesteps = linspace(flow.timespan(1),flow.timespan(2),...
    diff(flow.timespan)*framesPerSecond+1);

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
    delay = 1/framesPerSecond - toc;
    if delay < 0
        warning('animate_flow:negative_delay','Frame rate too high')
    end
    pause(delay)
    tic
end
