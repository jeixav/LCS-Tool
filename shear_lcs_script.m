function output = shear_lcs_script(input)

flow = set_flow_default(input.flow);

[flow,output.shearline] = compute_shear_lcs(flow,input.shearline);

output.flow = flow;
if isfield(input,'strainline')
    output.strainline = input.strainline;
end

mainAxes = setup_figure(flow);

showPlot = struct('shearlinePosFiltered',true,...
                  'shearlineNegFiltered',true);

plot_shear_lcs(mainAxes,flow,shearline,showPlot)

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
