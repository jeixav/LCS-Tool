function superminIndex = superminimize(position,value,superminDistance,...
    showPlot)
% SUPERMIN_LINE Superminimization on a line.

if nargin < 4
    showPlot = false;
end

% Find local minima by looking at nearest neighbours
localminIndex = false(size(position));
for iPosition = 2:(numel(position) - 1)
    if  (value(iPosition) < value(iPosition-1)) && ...
            (value(iPosition) < value(iPosition+1))
        localminIndex(iPosition) = true;
    end
end
localminIndex = find(localminIndex);

lmPosition = nan(size(position));
lmPosition(localminIndex) = position(localminIndex);
lmValue = nan(size(value));
lmValue(localminIndex) = value(localminIndex);

% Find superminima with a while loop as follows.
% 1. From all the local minima, find the global minimum and classify it as
% a superminimum.
% 2. Set all the local minima within ±superminimization_distance of this
% superminimum to NaN. These local minima will not be admissible as
% superminima.
% 3. Return to step 1 until all local minima have been set to NaN.
superminIndex = false(size(position));
while ~all(isnan(lmValue))
    [~,i] = min(lmValue);
    superminIndex(i) = true;
    min_interval = lmPosition(i) + .5*superminDistance*[-1 1];
    index = (lmPosition < min_interval(1)) | ...
        (lmPosition > min_interval(2));
    lmValue(~index) = nan;
end
superminIndex = find(superminIndex);

if showPlot
    if ~ishandle(showPlot)
        figure
        a1 = axes('nextplot','add','box','on','xgrid','on','ygrid','on');
    else
        a1 = showPlot;
    end
    plot(a1,position,value,'Marker','o','linestyle','none')
    plot(a1,position(superminIndex),value(superminIndex),'color','r',...
        'marker','o','linestyle','none','markerfacecolor','r')
end

end