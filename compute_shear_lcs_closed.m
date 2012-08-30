function [flow,shearline] = compute_shear_lcs_closed(flow,shearline)

if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    verbose.progressBar = true;
    verbose.stats = false;
    eigenvalueFromMainGrid = true;
   % [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
    %    eigenvalueFromMainGrid,verbose);
	
	flow.cgEigenvalue = zeros(flow.resolution(1) * flow.resolution(2), 2);
	flow.cgEigenvector = zeros(flow.resolution(1) * flow.resolution(2), 4);
end

%if ~all(isfield(shearline,{'etaPos','etaNeg','positionPos','positionNeg'}))
    verbose = true;
    %shearline = compute_shearline(flow,shearline,verbose);
	
	% convert tensor to array
	cgStrainTensor = zeros(flow.resolution(1) * flow.resolution(2), 4);
	%for n=1:length(flow.cgStrain)
	%	cgStrainTensor(n,1)=flow.cgStrain{n}(1,1);
	%	cgStrainTensor(n,2)=flow.cgStrain{n}(1,2);
	%	cgStrainTensor(n,3)=flow.cgStrain{n}(2,1);
	%	cgStrainTensor(n,4)=flow.cgStrain{n}(2,2);
	%end
	
	% call my function
	addpath('ClosedIntegralLines');
	flow.timespan
	flow.domain
	flow.resolution
	shearline.positionPos = []
	shearline.positionNeg = []
	[shearline.positionPos,shearline.positionNeg]=compute_closed_shearline(flow.timespan, flow.domain, flow.resolution, flow.cgEigenvalue, flow.cgEigenvector, cgStrainTensor);
	shearline.positionPos=shearline.positionPos';
	shearline.positionNeg=shearline.positionNeg';
	
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
