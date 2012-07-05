function input = ftle_script(input)
% Example:
% matlabpool('open')
% addpath('flow_templates')
% addpath('parfor_progress')
% travelingWave = traveling_wave;
% travelingWave = ftle_script(travelingWave);
% 
% To adjust the FTLE range:
% set(get(gca,'children'),'levelList',linspace(5e-4,.5,40))

input.flow = set_flow_default(input.flow);

method.name = 'eov';
method.params = struct('odeSolverOptions',input.flow.odeSolverOptions);
% method.name = 'fd';
% method.params = struct('eigenvalueFromMainGrid',true);
% input.flow.auxiliaryGridRelativeDelta = 5e-3;

verbose = struct('progressBar',true,'stats',true);

if ~all(isfield(input.flow,{'cgEigenvalue','cgEigenvector'}))
    [input.flow.cgEigenvalue input.flow.cgEigenvector] = ...
        eig_cgStrain(input.flow,method,verbose);
end

ftle = compute_ftle(input.flow.cgEigenvalue(:,2),...
    abs(diff(input.flow.timespan)));

axes = setup_figure(input.flow);
plot_ftle(axes,input.flow,ftle)

function a1 = setup_figure(flow)

figure
a1 = axes;
set(a1,'nextplot','add',...
    'box','on',...
    'DataAspectRatio',[1 1 1],...
    'DataAspectRatioMode','Manual',...
    'XGrid','on',...
    'YGrid','on',...
    'XLim',flow.domain(1,:),...
    'YLim',flow.domain(2,:))
xlabel(a1,'x')
ylabel(a1,'y')