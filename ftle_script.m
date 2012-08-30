% ftle_script Produce finite-time Lyapunov exponent plot
%
% DESCRIPTION
% input = ftle_script(input)
%
% EXAMPLE
% matlabpool('open')
% addpath('flow_templates')
% travelingWave = traveling_wave;
% travelingWave = ftle_script(travelingWave);
% 
% To adjust the FTLE range:
% set(get(gca,'children'),'levelList',linspace(5e-4,.5,40))

function input = ftle_script(input,verbose)

input.flow = set_flow_default(input.flow);

method.name = 'eov';
method.params = struct('odeSolverOptions',input.flow.odeSolverOptions);

if nargin == 1
    verbose = struct('progress',true,'stats',true);
end

if ~all(isfield(input.flow,{'cgEigenvalue','cgEigenvector'}))
    [input.flow.cgEigenvalue input.flow.cgEigenvector] = ...
        eig_cgStrain(input.flow,method,verbose);
end

ftle = compute_ftle(input.flow.cgEigenvalue(:,2),...
    abs(diff(input.flow.timespan)));

axes = setup_figure(input.flow);
plot_ftle(axes,input.flow,ftle)
