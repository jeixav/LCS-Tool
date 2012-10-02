function strainline = set_strainline_filtering_parameters(...
    filteringParameters,strainline)

strainline.filteringParameters = filteringParameters;

if isfield(strainline,'filteredSegmentIndex')
    strainline = rmfield(strainline,'filteredSegmentIndex');
end
