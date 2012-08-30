%% strain_lcs_script
% Compute and plot strainline LCSs
%
%% Description
% ...
%
%% Example
addpath('flow_templates')
matlabpool('open') % Optional; enable parallel computing.
doubleGyre = double_gyre;
doubleGyre = strain_lcs_script(doubleGyre);
