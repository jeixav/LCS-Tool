function [flow,shearline] = compute_shear_lcs_closed(flow,shearline)

if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    verbose.progress = true;
    verbose.stats = false;
    eigenvalueFromMainGrid = true;
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        eigenvalueFromMainGrid,verbose);
end

%if ~all(isfield(shearline,{'etaPos','etaNeg','positionPos','positionNeg'}))
    verbose = true;
    %shearline = compute_shearline(flow,shearline,verbose);
	
	addpath('ClosedIntegralLines');
	flow.timespan
	flow.domain
	flow.resolution
	shearline.positionPos = []
	shearline.positionNeg = []
	[shearline.positionPos,shearline.positionNeg]=compute_closed_shearline(flow.timespan, flow.domain, flow.resolution, flow.cgEigenvalue, flow.cgEigenvector);
	%shearline.positionNeg=shearline.positionPos
	shearline.positionPos=shearline.positionPos'
	shearline.positionNeg=shearline.positionNeg'
	
	shearline.geodesicDeviationPos=shearline.positionPos;
	shearline.geodesicDeviationNeg=shearline.positionNeg;
	shearline.averageGeodesicDeviationPos=rand(numel(shearline.positionPos),1);
	shearline.averageGeodesicDeviationNeg=rand(numel(shearline.positionNeg),1);
	shearline.filteredIndexPos=logical(ones(numel(shearline.positionPos),1));
	shearline.filteredIndexNeg=logical(ones(numel(shearline.positionNeg),1));
	shearline
%end

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
