function status = odeOutputFcn(time,y,flag,handles,frame_delay)

flow = getappdata(handles.LCSTool,'flow');

if strcmp(flag,'init')
    delete(get(handles.axesParticles,'children'))
    n = size(y,1)/2;
    if isfield(flow,'periodicBc')
        position = impose_periodic_bc([y(1:n) y(n+1:end)],...
            flow.domain,flow.periodicBc);
        y(1:n) = position(:,1);
        y(n+1:end) = position(:,2);
    end
        p1 = plot(handles.axesParticles,y(1:n),y(n+1:end));
    set(p1,'LineStyle','none',...
        'Marker','o',...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black',...
        'Tag','particles')
    setappdata(handles.axesParticles,'particles',p1)
    t1 = text(0,0,['t = ',num2str(time(1))]);
    set(t1,'Parent',handles.axesParticles,...
        'units','normalized',...
        'position',[.9 -.1 0],...
        'Tag','timeText')
    setappdata(handles.axesParticles,'timeText',t1)
    drawnow
    tic
    pause(frame_delay)
    set(handles.playButton,'UserData',1)
    set(handles.stopButton,'Enable','On')
    status = 0;
elseif strcmp(flag,'done')
    set(handles.stopButton,'Enable','Off')
    status = 0;
else
    if flow.animationTimespan
        timesteps = size(y,2);
        n = size(y,1)/2;
        for i = 1:timesteps
            if isfield(flow,'periodicBc')
                position = impose_periodic_bc([y(1:n,i) y(n+1:end,i)],...
                    flow.domain,flow.periodicBc);
                y(1:n,i) = position(:,1);
                y(n+1:end,i) = position(:,2);
            end
            set(findobj(handles.axesParticles,'Tag','particles'),...
                'xdata',y(1:n,i),'ydata',y(n+1:end,i))
            set(findobj(handles.axesParticles,'Tag','timeText'),...
               'string',['t = ',num2str(time(i),'%.1f')])
            drawnow
            % Pause when calculations complete faster than required
            pause(frame_delay - toc)
            tic
        end
    end
    if get(handles.playButton,'UserData') == 0
        status = 1;
    else
        status = 0;
    end
end
