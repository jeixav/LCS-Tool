function [input,hAxes] = shear_lcs_script(input,showPlot,verbose)

narginchk(1,3)

verboseDefault.progress = true;
verboseDefault.stats = true;
if nargin < 3
    verbose = [];
end
verbose = set_default(verbose,verboseDefault);

showPlotDefault.shearlinePosFiltered = true;
showPlotDefault.shearlineNegFiltered = true;
if nargin < 2
    showPlot = [];
end
showPlot = set_default(showPlot,showPlotDefault);

input.flow = set_flow_default(input.flow);

[input.flow,input.shearline] = compute_shear_lcs(input.flow,input.shearline,verbose);

hAxes = setup_figure(input.flow.domain);

plot_shear_lcs(hAxes,input.flow,input.shearline,showPlot)
