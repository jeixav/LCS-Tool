% strain_lcs_script Compute and plot strainline LCSs
%
% DESCRIPTION
% ...
%
% EXAMPLE
% addpath('flow_templates')
% matlabpool('open') % Optional; enable parallel computing.
% doubleGyre = double_gyre;
% doubleGyre = strain_lcs_script(doubleGyre);
%

function output = strain_lcs_script(input,showPlot,verbose)

narginchk(1,3)

if isa(input,'char')
    input = load_input_file(input);
end

output.flow = set_flow_default(input.flow);

if nargin < 3
    verbose.progress = true;
    verbose.stats = true;
    verbose.graphs = false;
end

if nargin < 2
    showPlot = [];
end
showPlot = set_showPlot_default(showPlot);

[output.flow,output.strainline] = compute_strain_lcs(output.flow,...
    input.strainline,verbose);

mainAxes = setup_figure(output.flow);

plot_strain_lcs(mainAxes,output.flow,output.strainline,showPlot)

if isfield(input,'shearline')
    output.shearline = input.shearline;
end

if isfield(input,'noStretchLine')
    output.noStretchLine= input.noStretchLine;
end

function showPlot = set_showPlot_default(showPlot)
% Default values for showPlot structure.

showPlotDefault = struct(...
    'quiver',false,...
    'geodesicDeviation',false,...
    'strainlineInitialCondition',false,...
    'strainline',false,...
    'strainlineSegment',false,...
    'ftle',false,...
    'superminLine',false,...
    'strainlineFiltered',true);

if ~isfield(showPlot,'quiver')
    showPlot.quiver = showPlotDefault.quiver;
end

if ~isfield(showPlot,'geodesicDeviation')
    showPlot.geodesicDeviation = showPlotDefault.geodesicDeviation;
end

if ~isfield(showPlot,'strainlineInitialCondition')
    showPlot.strainlineInitialCondition = ...
        showPlotDefault.strainlineInitialCondition;
end

if ~isfield(showPlot,'strainline')
    showPlot.strainline = showPlotDefault.strainline;
end

if ~isfield(showPlot,'strainlineSegment')
    showPlot.strainlineSegment = showPlotDefault.strainlineSegment;
end

if ~isfield(showPlot,'ftle')
    showPlot.ftle = showPlotDefault.ftle;
end

if ~isfield(showPlot,'superminLine')
    showPlot.superminLine = showPlotDefault.superminLine;
end

if ~isfield(showPlot,'strainlineFiltered')
    showPlot.strainlineFiltered = showPlotDefault.strainlineFiltered;
end
