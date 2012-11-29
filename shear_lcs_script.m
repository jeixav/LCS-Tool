function input = shear_lcs_script(input,showPlot,verbose)

narginchk(1,3)

verboseDefault = struct('progress',true,'stats',true);
showPlotDefault = struct('shearlinePosFiltered',true,...
    'shearlineNegFiltered',true);
    
if nargin < 3
    verbose = [];
end
verbose = set_default(verbose,verboseDefault);

if nargin < 2
    showPlot = [];
end
showPlot = set_default(showPlot,showPlotDefault);

input.flow = set_flow_default(input.flow);

[input.flow,input.shearline] = compute_shear_lcs(input.flow,...
    input.shearline,verbose);

mainAxes = setup_figure(input.flow.domain);

plot_shear_lcs(mainAxes,input.flow,input.shearline,showPlot)
