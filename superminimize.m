% SUPERMIN_LINE Superminimization on a line.
function superminIndex = superminimize(position,value,superminDistance,...
    showPlot)

% Flag to control if end points can be minima
allowEndpointMinima = true;

data = [position.' value.'];
validateattributes(data,{'double'},{'ncols',2})

if nargin < 4
    showPlot = false;
end

% Local minimization
% Check if a point is a minimum by comparing to its two nearest neighbors.
[sortedData, sortIndex] = sortrows(data);
nPoints = size(sortedData,1);
localminIndex = false(nPoints,1);
for iPosition = 2:(nPoints - 1)
    if  (sortedData(iPosition,2) <= sortedData(iPosition-1,2)) && ...
            (sortedData(iPosition,2) <= sortedData(iPosition+1,2))
        localminIndex(iPosition) = true;
    end
end
% End points require special treament since they have one nearest
% neighbour not two.
if size(sortedData,1) > 1
    if allowEndpointMinima
        if sortedData(1,2) <= sortedData(2,2)
            localminIndex(1) = true;
        end
        if sortedData(end,2) < sortedData(end-1,2)
            localminIndex(end) = true;
        end
    end
end
localminIndex = find(localminIndex);

% Second minimization
% Find superminima with a while loop as follows.
% 1. From all the local minima, find the global minimum and classify it as
% a superminimum.
% 2. Set all the local minima within ±superminimization_distance/2 of this
% superminimum to NaN. These local minima will not be admissible as
% superminima.
% 3. Return to step 1 until all local minima have been set to NaN.
lmPosition = nan(nPoints,1);
lmPosition(localminIndex) = sortedData(localminIndex,1);
lmValue = nan(size(value));
lmValue(localminIndex) = sortedData(localminIndex,2);
superminIndex = false(nPoints,1);
while ~all(isnan(lmValue))
    [~,i] = min(lmValue);
    superminIndex(i) = true;
    min_interval = lmPosition(i) + .5*superminDistance*[-1 1];
    index = (lmPosition < min_interval(1)) | ...
        (lmPosition > min_interval(2));
    lmValue(~index) = nan;
end
superminIndex = sortIndex(superminIndex);

if showPlot
    figure
    a1 = axes('nextplot','add','box','on','xgrid','on','ygrid','on');
    plot(a1,sortedData(:,1),sortedData(:,2),'Marker','o','linestyle','-')
    hLocalMin = plot(sortedData(localminIndex,1),sortedData(localminIndex,2));
    set(hLocalMin,'color','k')
    set(hLocalMin,'marker','o')
    set(hLocalMin,'markerfacecolor','k')
    set(hLocalMin,'linestyle','none')
    hSuperMin = plot(a1,position(superminIndex),value(superminIndex));
    set(hSuperMin,'color','r')
    set(hSuperMin,'marker','o')
    set(hSuperMin,'linestyle','none')
    set(hSuperMin,'markerfacecolor','r')
    legend([hLocalMin hSuperMin],'Local minimum','Super minimum')
    drawnow
end
