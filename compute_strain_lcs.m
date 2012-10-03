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
    verbose.progress = true;
    strainline = compute_strainline(flow,strainline,verbose);
end

if ~isfield(strainline,'geodesicDeviation')
    cgPosition = initialize_ic_grid(flow.resolution,flow.domain);
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
    cgPosition = initialize_ic_grid(flow.resolution,flow.domain);
    strainline.relativeStretching = relative_stretching(...
        strainline.position,strainline.segmentIndex,cgPosition,...
        flow.cgEigenvalue(:,1),flow.resolution);
end

if ~isfield(strainline,'filteredSegmentIndex')
    switch strainline.filteringMethod
        case 'superminimization'
            plotSuperminLine = false;
            strainline.filteredSegmentIndex = superminimize_grid(...
                strainline.position,strainline.segmentIndex,...
                strainline.relativeStretching,...
                strainline.filteringParameters.distance,...
                flow.domain,strainline.filteringParameters.resolution,...
                plotSuperminLine);
        case 'hausdorff'
                strainline = hausdorff_filtering(strainline);
        otherwise
            error('Filtering method invalid')
    end
end