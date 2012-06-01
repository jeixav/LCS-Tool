function output = plot_hyperbolic_lcs(input,showPlot)

if isa(input,'char')
    input = load_input_file(input);
end

flow = set_flow_default(input.flow);
strainline = input.strainline;

if nargin == 1
    showPlot = set_showPlot_default([]);
elseif nargin == 2
    showPlot = set_showPlot_default(showPlot);
else
    error('Incorrect number of input arguments.')
end

if ~isfield(flow,'finalPosition')
    initialPosition = initial_position(flow.domain,flow.resolution);
    deltaX = (flow.domain(1,2) - flow.domain(1,1))...
        /double(flow.resolution(1))*flow.auxiliaryGridRelativeDelta;
    delta = deltaX;
    auxiliaryPosition = auxiliary_position(initialPosition,delta);
    
    flow.odeSolverOptions = odeset(flow.odeSolverOptions,'outputFcn',...
        @ode_progress_bar);
    flow.finalPosition = integrate_flow(flow,auxiliaryPosition);
end

if ~isfield(flow,'cgStrain')
    flow.cgStrain = compute_cgStrain(flow.finalPosition,...
        delta);
end

if ~all(isfield(flow,{'cgEigenvector','cgEigenvector'}))
    [flow.cgEigenvector,flow.cgEigenvalue] = ...
        arrayfun(@eig_array,flow.cgStrain(:,1),...
        flow.cgStrain(:,2),flow.cgStrain(:,3),...
        'UniformOutput',false);
    flow.cgEigenvalue = cell2mat(flow.cgEigenvalue);
    flow.cgEigenvector = cell2mat(flow.cgEigenvector);
end

if ~isfield(strainline,'position')
    scaleXi1 = false;
    cgPosition = initial_position(flow.domain,flow.resolution);
    strainline.odeSolverOptions = odeset('outputFcn',@ode_progress_bar);
    strainline.position = compute_strainline(flow,...
        strainline,cgPosition,flow.cgEigenvalue,...
        flow.cgEigenvector(:,1:2),scaleXi1);
end

if ~isfield(strainline,'geodesicDeviation')
    cgPosition = initial_position(flow.domain,flow.resolution);
    strainline.geodesicDeviation = geodesic_deviation_strainline(...
        strainline.position,cgPosition,flow.cgEigenvalue(:,2),...
        flow.cgEigenvector,flow.resolution);
end
geodesic_deviation_stats(strainline.geodesicDeviation,true);

% Graphics
mainAxes = setup_figure(flow);

plot_ftle(mainAxes,flow)
if ~isfield(showPlot,'ftle') || showPlot.ftle == false
    set(findobj(mainAxes,'tag','ftle'),'visible','off')
end

cgPosition = initial_position(flow.domain,flow.resolution);
quiver(mainAxes,cgPosition(:,1),cgPosition(:,2),...
    flow.cgEigenvector(:,1),flow.cgEigenvector(:,2),'tag','quiver')
if ~isfield(showPlot,'quiver') || showPlot.quiver == false;
    set(findobj(mainAxes,'tag','quiver'),'visible','off')
end

plot_strainline(mainAxes,strainline)
if ~isfield(showPlot,'strainline') || showPlot.strainline == false
    set(findobj(mainAxes,'tag','strainline'),'visible','off')
end

strainlineIc = initialize_ic_grid(strainline.resolution,flow.domain);
plot(mainAxes,strainlineIc(:,1),strainlineIc(:,2),...
    'MarkerFaceColor','k',...
    'MarkerEdgeColor','k',...
    'Marker','o',...
    'LineStyle','none',...
    'Tag','strainlineInitialCondition')
if ~isfield(showPlot,'strainlineInitialCondition') || ...
        showPlot.strainlineInitialCondition == false
    set(findobj(mainAxes,'tag','strainlineInitialCondition'),'visible',...
        'off')
end

if showPlot.geodesicDeviation
    geodesicDeviationIndex = cellfun(...
        @(input)input<strainline.geodesicDeviationTol,...
        strainline.geodesicDeviation,'UniformOutput',false);
    cellfun(@(position,index)plot_geodesic_deviation_points(...
        position,index,mainAxes),strainline.position,geodesicDeviationIndex)
end
    
if ~isfield(strainline,'segmentIndex')
    strainline.segmentIndex = find_segments(strainline.position,...
        strainline.geodesicDeviation,...
        strainline.geodesicDeviationTol,...
        strainline.lengthTol);
    nSegments = sum(cellfun(@(input)size(input,1),strainline.segmentIndex));
    disp(['Number of strainline segments: ',num2str(nSegments)])
end

plot_strainline_segment(mainAxes,strainline.position,...
    strainline.segmentIndex)
if ~isfield(showPlot,'strainlineSegment') || ...
        showPlot.strainlineSegment == false
    set(findobj(mainAxes,'tag','strainlineSegment'),'visible','off')
end

if ~isfield(strainline,'relativeStretching')
    cgPosition = initial_position(flow.domain,flow.resolution);
    strainline.relativeStretching = relative_stretching(...
        strainline.position,strainline.segmentIndex,cgPosition,...
        flow.cgEigenvalue(:,1),flow.resolution);
end

if showPlot.strainlineFiltered || showPlot.strainlineSegment
    
    switch strainline.filteringMethod
        case 'hausdorff'

            if showPlot.superminLine == true
                showPlot.superminLine = mainAxes;
            end
            
            if ~all(isfield(strainline,{'filteredStrainlineIndex',...
                    'hausdorffDistance'}))
                strainline = hausdorff_filtering(strainline,...
                    showPlot.superminLine);
            end
            nSegments = sum(cellfun(@sum,strainline.filteredSegmentIndex));
            fprintf('Number of LCS segments: %g\n',nSegments)

        case 'superminimization'

            if showPlot.superminLine == true
                showPlot.superminLine = mainAxes;
            end
            
            if ~isfield(strainline,'filteredStrainlineIndex')
                strainline.filteredStrainlineIndex = superminimize_grid(...
                    strainline.position,strainline.segmentIndex,...
                    strainline.relativeStretching,...
                    strainline.superminimizationDistance,...
                    flow.domain,strainline.resolution,showPlot.superminLine);
            end
            % fprintf(['Done (',num2str(toc(filteringStrainlineTime)),'s.)\n'])
    
        otherwise
            error('Filtering method not recognized')
    end

    if showPlot.strainlineFiltered
        plot_filtered_strainline(mainAxes,strainline.position,...
            strainline.segmentIndex,strainline.filteredSegmentIndex)
    end
    
end

output.flow = flow;
output.strainline = strainline;
if isfield(input,'shearline')
    output.shearline = input.shearline;
end

function a1 = setup_figure(flow)

figure
a1 = axes;
set(a1,'nextplot','add',...
    'box','on',...
    'DataAspectRatio',[1 1 1],...
    'DataAspectRatioMode','Manual',...
    'XGrid','on',...
    'YGrid','on',...
    'XLim',flow.domain(1,:),...
    'YLim',flow.domain(2,:))
xlabel(a1,'x')
ylabel(a1,'y')

function showPlot = set_showPlot_default(showPlot)
% Add default values to showPlot structure.

showPlotDefault = struct(...
    'quiver',false,...
    'geodesicDeviation',false,...
    'strainlineInitialCondition',false,...
    'strainline',false,...
    'strainlineSegment',false,...
    'ftle',false,...
    'superminLine',false,...
    'strainlineFiltered',true);

if ~isfield(showPlot,'quiver')
    showPlot.quiver = showPlotDefault.quiver;
end

if ~isfield(showPlot,'geodesicDeviation')
    showPlot.geodesicDeviation = showPlotDefault.geodesicDeviation;
end

if ~isfield(showPlot,'strainlineInitialCondition')
    showPlot.strainlineInitialCondition = ...
        showPlotDefault.strainlineInitialCondition;
end

if ~isfield(showPlot,'strainline')
    showPlot.strainline = showPlotDefault.strainline;
end

if ~isfield(showPlot,'strainlineSegment')
    showPlot.strainlineSegment = showPlotDefault.strainlineSegment;
end

if ~isfield(showPlot,'ftle')
    showPlot.ftle = showPlotDefault.ftle;
end

if ~isfield(showPlot,'superminLine')
    showPlot.superminLine = showPlotDefault.superminLine;
end

if ~isfield(showPlot,'strainlineFiltered')
    showPlot.strainlineFiltered = showPlotDefault.strainlineFiltered;
end

function plot_geodesic_deviation_points(position,index,axes)
plot(axes,position(index,1),position(index,2),...
    'MarkerEdgeColor','r',...
    'Marker','o',...
    'LineStyle','none',...
    'Tag','geodesicDeviationPoint')

function plot_strainline(axes,strainline)

cellfun(@(position)plot(axes,position(:,1),position(:,2),'color','k',...
    'Tag','strainline'),strainline.position)

function plot_ftle(axes,flow)

ftle = compute_ftle(flow.cgEigenvalue(:,2),abs(diff(flow.timespan)));
ftle = reshape(ftle,fliplr(flow.resolution));
h = pcolor(axes,...
    linspace(flow.domain(1,1),flow.domain(1,2),flow.resolution(1)),...
    linspace(flow.domain(2,1),flow.domain(2,2),flow.resolution(2)),...
    ftle);
set(h,'tag','ftle')
shading(axes,'interp')
