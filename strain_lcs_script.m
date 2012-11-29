% strain_lcs_script Compute and plot strainline LCSs
%
% EXAMPLE
% addpath('flow_templates')
% matlabpool('open') % Optional; enable parallel computing.
% doubleGyre = double_gyre;
% doubleGyre = strain_lcs_script(doubleGyre);

function output = strain_lcs_script(input,showPlot,verbose)

% FIXME Code block identical with strain_lcs_script
narginchk(1,3)

verboseDefault = struct('graphs',true,'progress',true,'stats',true);
if nargin < 3
    verbose = [];
end
verbose = set_default(verbose,verboseDefault);

showPlotDefault = struct(...
    'quiver',false,...
    'geodesicDeviation',false,...
    'strainlineInitialCondition',false,...
    'strainline',false,...
    'strainlineSegment',false,...
    'ftle',false,...
    'superminLine',false,...
    'strainlineFiltered',true);
if nargin < 2
    showPlot = [];
end
showPlot = set_default(showPlot,showPlotDefault);

output.flow = input.flow;

[output.flow,output.strainline] = compute_strain_lcs(output.flow,...
    input.strainline,verbose);

mainAxes = setup_figure(output.flow.domain);

plot_strain_lcs(mainAxes,output.flow,output.strainline,showPlot)

if isfield(input,'shearline')
    output.shearline = input.shearline;
end

if isfield(input,'noStretchLine')
    output.noStretchLine= input.noStretchLine;
end
