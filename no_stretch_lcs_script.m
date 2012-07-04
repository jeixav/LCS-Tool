function input = no_stretch_lcs_script(input)

input.flow = set_flow_default(input.flow);

[input.flow,input.noStretchLine] = compute_no_stretch_lcs(input.flow,...
    input.noStretchLine);

mainAxes = setup_figure(input.flow);

showPlot = struct('noStretchLinePos',true,'noStretchLineNeg',true);

plot_no_stretch_lcs(mainAxes,input.flow,input.noStretchLine,showPlot)

function a1 = setup_figure(flow)

figure
a1 = axes;
set(a1,'nextplot','add',...
    'box','on',...
    'DataAspectRatio',[1 1 1],...
    'DataAspectRatioMode','Manual',...
    'XGrid','on',...
    'YGrid','on',...
    'XLim',flow.domain(1,:),...
    'YLim',flow.domain(2,:))
xlabel(a1,'x')
ylabel(a1,'y')
