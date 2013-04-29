function hPlot = plot_filtered_strainline(axes,position,segmentIndex,superminIndex)

nStrainlines = length(superminIndex);

hPlot = nan(nStrainlines,1);

for iStrainline = 1:nStrainlines
    superminIndexLocal = find(superminIndex{iStrainline});
    if superminIndexLocal
        xLocal = position{iStrainline}(:,1);
        yLocal = position{iStrainline}(:,2);
        startIndex = segmentIndex{iStrainline}(superminIndexLocal,1);
        stopIndex = segmentIndex{iStrainline}(superminIndexLocal,2);
        hPlot(iStrainline) = arrayfun(@(start,stop)plotArray(axes,xLocal,yLocal,start,stop),startIndex,stopIndex);
    end
end

hPlot = hPlot(~isnan(hPlot));

function hPlot = plotArray(axes,xLocal,yLocal,startIndex,stopIndex)

hPlot = plot(axes,xLocal(startIndex:stopIndex),yLocal(startIndex:stopIndex));
set(hPlot,'Color','r')
set(hPlot,'Tag','strainlineFiltered')
