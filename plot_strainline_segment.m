function plot_strainline_segment(axes,strainlinePosition,segmentIndex)

nStrainlines = size(segmentIndex,2);
for iStrainline = 1:nStrainlines
    nSegments = size(segmentIndex{iStrainline},1);
    for iSegment = 1:nSegments
        startIndex = segmentIndex{iStrainline}(iSegment,1);
        stopIndex = segmentIndex{iStrainline}(iSegment,2);
        strainlinePositionX = strainlinePosition{iStrainline}...
            (startIndex:stopIndex,1);
        strainlinePositionY = strainlinePosition{iStrainline}...
            (startIndex:stopIndex,2);
        plot(axes,strainlinePositionX,strainlinePositionY,...
            'color','r','tag','strainlineSegment')
    end
end

end
