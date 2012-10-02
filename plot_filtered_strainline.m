function plot_filtered_strainline(axes,position,segmentIndex,superminIndex)

nStrainlines = length(superminIndex);

for iStrainline = 1:nStrainlines
    superminIndexLocal = find(superminIndex{iStrainline});
    if superminIndexLocal
        xLocal = position{iStrainline}(:,1);
        yLocal = position{iStrainline}(:,2);
        startIndex = segmentIndex{iStrainline}(superminIndexLocal,1);
        stopIndex = segmentIndex{iStrainline}(superminIndexLocal,2);
        arrayfun(@(start,stop) plotArray(axes,xLocal,yLocal,start,stop),...
            startIndex,stopIndex)
    end
end

end

function plotArray(axes,xLocal,yLocal,startIndex,stopIndex)

plot(axes,xLocal(startIndex:stopIndex),yLocal(startIndex:stopIndex),...
    'Color','r','Tag','filteredStrainline')

end