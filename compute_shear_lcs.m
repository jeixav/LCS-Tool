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
        warning([mfilename,':defaultcgStrainMethodName'],...
            ['flow.cgStrainMethod.name not set; using default: ',...
            cgStrainMethod.name])
    else
        cgStrainMethod = flow.cgStrainMethod;
    end
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        cgStrainMethod,verbose);
end

if ~all(isfield(shearline,{'etaPos','etaNeg','positionPos','positionNeg'}))
    shearline = compute_shearline(flow,shearline,verbose);
end

if ~all(isfield(shearline,{'geodesicDeviationPos','geodesicDeviationNeg',...
        'averageGeodesicDeviationPos','averageGeodesicDeviationNeg'}))
    shearline = geodesic_deviation_shearline(flow,shearline,...
        verbose.progress);
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
