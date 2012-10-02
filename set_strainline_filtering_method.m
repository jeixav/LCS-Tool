function strainline = set_strainline_filtering_method(...
    filteringMethod,strainline)

strainline.filteringMethod = filteringMethod;

if isfield(strainline,'filteredSegmentIndex')
    strainline = rmfield(strainline,'filteredSegmentIndex');
end
