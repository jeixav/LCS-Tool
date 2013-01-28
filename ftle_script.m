% ftle_script Produce finite-time Lyapunov exponent plot
%
% DESCRIPTION
% input = ftle_script(input)
%
% EXAMPLE
% matlabpool('open')
% pctRunOnAll javaaddpath('ParforProgress2')
% addpath('flow_templates')
% travelingWave = traveling_wave;
% travelingWave.flow = ftle_script(travelingWave.flow);
% travelingWave = strain_lcs_script(travelingWave);
% close
% plot_filtered_strainline(gca,travelingWave.strainline.position,...
%     travelingWave.strainline.segmentIndex,...
%     travelingWave.strainline.filteredSegmentIndex)
% set(findobj(gca,'tag','strainlineFiltered'),'linewidth',1)
% 
% To adjust the FTLE range:
% set(get(gca,'children'),'levelList',linspace(5e-4,.5,40))

function flow = ftle_script(flow,verbose)

% input.flow = set_flow_default(input.flow);

% method.name = 'eov';
% method.params = struct('odeSolverOptions',input.flow.odeSolverOptions);

if nargin == 1
    verbose = struct('progress',true,'stats',true);
end

% FIXME This if-statement is identical with one in compute_strain_lcs
if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    if ~isfield(flow,'cgStrainMethod')
        cgStrainMethod.name = 'equationOfVariation';
        warning([mfilename,':defaultcgStrainMethodName'],...
            ['flow.cgStrainMethod.name not set; using default: ',...
            cgStrainMethod.name])
    else
        cgStrainMethod = flow.cgStrainMethod;
    end
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        cgStrainMethod,[],verbose);
end

ftle = compute_ftle(flow.cgEigenvalue(:,2),abs(diff(flow.timespan)));

axes = setup_figure(flow.domain);
plot_ftle(axes,flow,ftle)
