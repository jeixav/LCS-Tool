function input = shear_lcs_script(input)

input.flow = set_flow_default(input.flow);

[input.flow,input.shearline] = compute_shear_lcs(input.flow,...
    input.shearline);

mainAxes = setup_figure(input.flow);

showPlot = struct('shearlinePosFiltered',true,...
                  'shearlineNegFiltered',true,...
                  'shearlinePosClosed',true,...
                  'shearlineNegClosed',true);

plot_shear_lcs(mainAxes,input.flow,input.shearline,showPlot)
