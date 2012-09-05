% animate_flow Display flow animation
%
% DESCRIPTION
% flow = animate_flow(flow,framesPerSecond,verbose)
%
% EXAMPLE
% addpath('flow_templates')
% doubleGyre = double_gyre;
% doubleGyre.flow = animate_flow(doubleGyre.flow);

function flow = animate_flow(flow,framesPerSecond,verbose)

if nargin < 3
    verbose.progress = true;
end

if nargin < 2
    framesPerSecond = 10;
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
    if isfield(flow,'symDerivative') && ~isfield(flow,'derivative')
        flow.derivative = sym2fun(flow.symDerivative);
    end
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
