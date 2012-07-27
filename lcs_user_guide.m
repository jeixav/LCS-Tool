%% Lagrangian Coherent Structures Toolbox User Guide
%

%% Flow definition
% The minimum required to define a flow is...

%% Example 1: Flow animation
%
addpath('flow_templates')
doubleGyre = double_gyre;
doubleGyre.flow = animate_flow(doubleGyre.flow);

%% Example 2: Hyperbolic barriers
%
addpath('flow_templates')
matlabpool('open')  % Optional; enable parallel computing
doubleGyre = double_gyre;
doubleGyre = strain_lcs_script(doubleGyre);

%% Example 3: Parabolic barriers
% 

%% Example 4: Elliptic barriers
%
