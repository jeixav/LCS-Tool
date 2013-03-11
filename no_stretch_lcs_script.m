function input = no_stretch_lcs_script(input)

[input.flow,input.noStretchLine] = compute_no_stretch_lcs(input.flow,...
    input.noStretchLine);

mainAxes = setup_figure(input.flow.domain);

showPlot = struct('noStretchLinePos',true,'noStretchLineNeg',true);

plot_no_stretch_lcs(mainAxes,input.flow,input.noStretchLine,showPlot)
