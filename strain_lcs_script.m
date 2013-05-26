% strain_lcs_script Compute and plot strainline LCSs
%
% SYNTAX
% [input,hAxes] = strain_lcs_script(input,showPlot,verbose)
%
% EXAMPLE
% addpath('flow_templates')
% matlabpool('open') % Optional; enable parallel computing.
% doubleGyre = double_gyre;
% doubleGyre = strain_lcs_script(doubleGyre);

function [input,hAxes] = strain_lcs_script(input,showPlot,verbose)

narginchk(1,3)

verboseDefault.graphs = true;
verboseDefault.progress = false;
verboseDefault.stats = true;
if nargin < 3
    verbose = [];
end
verbose = set_default(verbose,verboseDefault);

showPlotDefault.quiver = false;
showPlotDefault.geodesicDeviation = false;
showPlotDefault.strainlineInitialCondition = false;
showPlotDefault.strainline = false;
showPlotDefault.strainlineSegment = false;
showPlotDefault.ftle = false;
showPlotDefault.superminLine = false;
showPlotDefault.strainlineFiltered = true;
if nargin < 2
    showPlot = [];
end
showPlot = set_default(showPlot,showPlotDefault);

[input.flow,input.strainline] = compute_strain_lcs(input.flow,input.strainline,verbose);

hAxes = setup_figure(input.flow.domain);

plot_strain_lcs(hAxes,input.flow,input.strainline,showPlot)
