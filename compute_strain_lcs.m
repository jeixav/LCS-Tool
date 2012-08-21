function [flow,strainline] = compute_strain_lcs(flow,strainline,verbose)

if nargin < 3
    verbose.progress = false;
end

if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    verbose.progress = true;
    verbose.stats = false;
    cgStrainMethod.name = 'eov';
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        cgStrainMethod,verbose);
end

if ~isfield(strainline,'position')
    strainline = compute_strainline(flow,strainline);
end

if ~isfield(strainline,'geodesicDeviation')
    cgPosition = initial_position(flow.domain,flow.resolution);
    strainline.geodesicDeviation = geodesic_deviation_strainline(...
        strainline.position,cgPosition,flow.cgEigenvalue(:,2),...
        flow.cgEigenvector,flow.resolution);
end

geodesic_deviation_stats(strainline.geodesicDeviation,true);

if ~isfield(strainline,'segmentIndex')
    strainline.segmentIndex = find_segments(strainline.position,...
        strainline.geodesicDeviation,...
        strainline.geodesicDeviationTol,...
        strainline.lengthTol);
    nSegments = sum(cellfun(@(input)size(input,1),strainline.segmentIndex));
    disp(['Number of strainline segments: ',num2str(nSegments)])
end

if ~isfield(strainline,'relativeStretching')
    cgPosition = initial_position(flow.domain,flow.resolution);
    strainline.relativeStretching = relative_stretching(...
        strainline.position,strainline.segmentIndex,cgPosition,...
        flow.cgEigenvalue(:,1),flow.resolution);
end

if ~all(isfield(strainline,{'hausdorffDistance','filteredSegmentIndex'}))
    strainline = hausdorff_filtering(strainline);
    nSegments = sum(cellfun(@sum,strainline.filteredSegmentIndex));
    fprintf('Number of LCS segments: %g\n',nSegments)
end
