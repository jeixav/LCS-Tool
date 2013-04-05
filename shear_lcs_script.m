% shear_lcs_script Compute and plot shear LCSs
%
% SYNTAX
% [input,hAxes] = shear_lcs_script(input)
% [input,hAxes] = shear_lcs_script(input,showPlot)
% [input,hAxes] = shear_lcs_script(input,showPlot,verbose)

function [input,hAxes] = shear_lcs_script(input,varargin)

narginchk(1,3)

p = inputParser;
p.StructExpand = false;

validationFcn = @(x) isstruct(x) && isfield(x,'flow') && isfield(x,'shearline');
addRequired(p,'input',validationFcn)

showPlotDefault.shearlinePosFiltered = true;
showPlotDefault.shearlineNegFiltered = true;
addOptional(p,'showPlot',showPlotDefault,@isstruct)

verboseDefault.progress = true;
verboseDefault.stats = true;
addOptional(p,'verbose',verboseDefault,@isstruct)

parse(p,input,varargin{:})

input = p.Results.input;
showPlot = p.Results.showPlot;
verbose = p.Results.verbose;

[input.flow,input.shearline] = compute_shear_lcs(input.flow,input.shearline,verbose);

hAxes = setup_figure(input.flow.domain);

plot_shear_lcs(hAxes,input.flow,input.shearline,showPlot)
