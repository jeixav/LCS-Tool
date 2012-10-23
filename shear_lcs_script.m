function input = shear_lcs_script(input,showPlot)

if nargin == 1
    showPlot = struct('shearlinePosFiltered',true,...
        'shearlineNegFiltered',true);
end

input.flow = set_flow_default(input.flow);

[input.flow,input.shearline] = compute_shear_lcs(input.flow,...
    input.shearline);

mainAxes = setup_figure(input.flow);

plot_shear_lcs(mainAxes,input.flow,input.shearline,showPlot)
