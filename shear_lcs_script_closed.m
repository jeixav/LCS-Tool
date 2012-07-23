function output = shear_lcs_script_closed(input)

flow = set_flow_default(input.flow);

[output.flow,output.shearline] = compute_shear_lcs_closed(flow,input.shearline);

if isfield(input,'strainline')
    output.strainline = input.strainline;
end

mainAxes = setup_figure(output.flow);

showPlot = struct('shearlinePosFiltered',true,...
                  'shearlineNegFiltered',true);

plot_shear_lcs(mainAxes,output.flow,output.shearline,showPlot)
