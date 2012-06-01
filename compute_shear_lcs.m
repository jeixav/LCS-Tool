function [flow,shearline] = compute_shear_lcs(flow,shearline)

if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    verbose.progressBar = true;
    verbose.stats = false;
    eigenvalueFromMainGrid = true;
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        eigenvalueFromMainGrid,verbose);
end

if ~all(isfield(shearline,{'etaPos','etaNeg','positionPos','positionNeg'}))
    verbose = true;
    shearline = compute_shearline(flow,shearline,verbose);
end

if ~all(isfield(shearline,{'geodesicDeviationPos','geodesicDeviationNeg',...
        'averageGeodesicDeviationPos','averageGeodesicDeviationNeg'}))
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

index = averageGeodesicDeviation <= tol;
