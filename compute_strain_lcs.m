function [flow,strainline] = compute_strain_lcs(flow,strainline,verbose)

if nargin < 3
    verbose.progress = true;
    verbose.stats = true;
    verbose.graphs = false;
end

if ~all(isfield(flow,{'cgEigenvalue','cgEigenvector'}))
    if ~isfield(flow,'cgStrainMethod')
        cgStrainMethod.name = 'equationOfVariation';
        warning('compute_strain_lcs:defaultcgStrainMethodName',...
            ['flow.cgStrainMethod.name not set; using default: ',...
            cgStrainMethod.name])
    else
        cgStrainMethod = flow.cgStrainMethod;
    end
    [flow.cgEigenvalue,flow.cgEigenvector] = eig_cgStrain(flow,...
        cgStrainMethod,verbose);
end

if ~isfield(strainline,'position')
    strainline = compute_strainline(flow,strainline,verbose);
end

if ~isfield(strainline,'geodesicDeviation')
    cgPosition = initialize_ic_grid(flow.resolution,flow.domain);
    strainline.geodesicDeviation = geodesic_deviation_strainline(...
        strainline.position,cgPosition,flow.cgEigenvalue(:,2),...
        flow.cgEigenvector,flow.resolution);
    geodesic_deviation_stats(strainline.geodesicDeviation,true);
end

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
        case 'superminimize'
            plotSuperminLine = verbose.graphs;
            matlabpoolClosed = false;
            if plotSuperminLine && matlabpool('size')
                warning('compute_strain_lcs:plotSuperminLine',...
                    ['plotSuperminLine does not work when matlabpool '...
                    'is in use. Temporarily closing matlabpool.'])
                matlabpool('close')
                matlabpoolClosed = true;
            end
            strainline.filteredSegmentIndex = superminimize_grid(...
                strainline.position,strainline.segmentIndex,...
                strainline.relativeStretching,...
                strainline.filteringParameters.distance,...
                flow.domain,strainline.filteringParameters.resolution,...
                plotSuperminLine);
            if matlabpoolClosed
                matlabpool('open')
            end
        case 'hausdorff'
            strainline = hausdorff_filtering(strainline);
        case 'minimin'
            strainline = minimin_filtering(strainline);
        otherwise
            error('strainline.filteringMethod invalid')
    end
    nfilteredSegment = sum(cellfun(@sum,strainline.filteredSegmentIndex));
    disp(['Number of filtered segments: ' num2str(nfilteredSegment)])
end
