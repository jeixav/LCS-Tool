% ftle_script Produce finite-time Lyapunov exponent plot
%
% DESCRIPTION
% input = ftle_script(input)
%
% EXAMPLE
% pctRunOnAll javaaddpath('ParforProgress2')
% addpath('flow_templates')
% travelingWave = traveling_wave;
% travelingWave.flow = ftle_script(travelingWave.flow);
% travelingWave = strain_lcs_script(travelingWave);
% close
% plot_filtered_strainline(gca,travelingWave.strainline.position,...
%     travelingWave.strainline.segmentIndex,...
%     travelingWave.strainline.filteredSegmentIndex);
% set(findobj(gca,'tag','strainlineFiltered'),'linewidth',1)
% 
% To adjust the FTLE range:
% caxis(5e-4,.5)

function flow = ftle_script(flow,verbose)

if nargin == 1
    verbose = true;
end

% FIXME This if-statement is identical with one in compute_strain_lcs and
% compute_shear_lcs
if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    if ~isfield(flow,'cgStrainMethod')
        cgStrainMethod.name = 'equationOfVariation';
        warning([mfilename,':defaultCgStrainMethodName'],['flow.cgStrainMethod.name not set; using default: ',cgStrainMethod.name])
    else
        cgStrainMethod = flow.cgStrainMethod;
    end
    
    if ~isfield(flow,'cgStrainCustomEigMethod')
        cgStrainCustomEigMethod = false;
        warning([mfilename,':defaultCgStrainCustomEigMethod'],['flow.cgStrainCustomEigMethod not set; using default: ',num2str(cgStrainCustomEigMethod)])
    else
        cgStrainCustomEigMethod = flow.cgStrainCustomEigMethod;
    end
    
    if ~isfield(flow,'coupledIntegration')
        coupledIntegration = false;
    else
        coupledIntegration = flow.coupledIntegration;
    end
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,cgStrainMethod,cgStrainCustomEigMethod,coupledIntegration,verbose);
end

ftle = compute_ftle(flow.cgEigenvalue(:,2),abs(diff(flow.timespan)));

axes = setup_figure(flow.domain);
plot_ftle(axes,flow,ftle);
