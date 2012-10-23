function superminIndex = superminimize2(position,value,superminDistance,...
    showPlot)

if nargin < 4
    showPlot = false;
end

if showPlot
    % Save original value
    valueOrig = value;
end

superminIndex = nan(size(position));
counterIndex = 0;

while ~all(isnan(value))
    counterIndex = counterIndex + 1;
    [~,minIndex] = min(value);
    value(minIndex) = nan;
    minInterval = position(minIndex) + .5*superminDistance*[-1 1];

    % CleanSuperminIndex is the same as superminIndex, but with NaN values
    % removed.
    cleanSuperminIndex = superminIndex(~isnan(superminIndex));
    if ~any((position(cleanSuperminIndex) > minInterval(1)) & ...
            (position(cleanSuperminIndex) <= minInterval(2)))
        superminIndex(counterIndex) = minIndex;
    end
end

superminIndex = superminIndex(~isnan(superminIndex));

if showPlot
    if ~ishandle(showPlot)
        figure
        a1 = axes('nextplot','add','box','on','xgrid','on','ygrid','on');
    else
        a1 = showPlot;
    end
    sortedData = sortrows([position.' ,valueOrig.']);
    plot(a1,sortedData(:,1),sortedData(:,2),'Marker','o','linestyle','-');
    plot(a1,position(superminIndex),valueOrig(superminIndex),...
        'color','r','marker','o','linestyle','none','markerfacecolor','r');
    drawnow
end
