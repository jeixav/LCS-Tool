%% animate_flow
% Display flow animation
%
%% Description
% flow = animate_flow(flow,framesPerSecond,verbose)
%
%% Example
addpath('flow_templates')
doubleGyre = double_gyre;
doubleGyre.flow = animate_flow(doubleGyre.flow);