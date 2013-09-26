function [flow,shearline] = compute_shear_lcs(flow,shearline,verbose)

narginchk(2,3)

verboseDefault = struct('progress',true,'stats',true);
if nargin < 3
    verbose = [];
end
verbose = set_default(verbose,verboseDefault);

% FIXME This if-statement is identical with one in compute_strain_lcs
if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))

    if ~isfield(flow,'cgStrainMethod')
        cgStrainMethod.name = 'equationOfVariation';
    else
        cgStrainMethod = flow.cgStrainMethod;
    end
    
    if ~isfield(flow,'cgStrainCustomEigMethod')
        cgStrainCustomEigMethod = false;
    else
        cgStrainCustomEigMethod = flow.cgStrainCustomEigMethod;
    end
    
    if ~isfield(flow,'coupledIntegration')
        coupledIntegration = false;
    else
        coupledIntegration = flow.coupledIntegration;
    end

    [flow.cgEigenvalue,flow.cgEigenvector,flow.cgStrain] = eig_cgStrain(flow,cgStrainMethod,cgStrainCustomEigMethod,coupledIntegration,verbose.stats);
end

if ~all(isfield(shearline,{'etaPos','etaNeg','positionPos','positionNeg'}))
    shearline = compute_shearline(flow,shearline,verbose);
end

if ~all(isfield(shearline,{'geodesicDeviationPos','geodesicDeviationNeg','averageGeodesicDeviationPos','averageGeodesicDeviationNeg'}))
    shearline = geodesic_deviation_shearline(flow,shearline);
end

if ~all(isfield(shearline,{'filteredIndexPos','filteredIndexNeg'}))
    shearline.filteredIndexPos = filter_shearline(...
        shearline.averageGeodesicDeviationPos,...
        shearline.averageGeodesicDeviationPosTol);
    shearline.filteredIndexNeg = filter_shearline(...
        shearline.averageGeodesicDeviationNeg,...
        shearline.averageGeodesicDeviationNegTol);
end

function index = filter_shearline(averageGeodesicDeviation,tol)

% FIXME Workaround while boundary element geodesic is not calculated and
% yields NaNs
if isinf(tol)
    warning('averageGeodesicDeviation = inf workaround used')
    index = true(size(averageGeodesicDeviation));
else
    index = averageGeodesicDeviation <= tol;
end
