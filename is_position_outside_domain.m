function output = is_position_outside_domain(position,domain)

nStrainlines = size(position,2);
positionXMin = min(arrayfun(@(i)min(position{i}(:,1)),1:nStrainlines));
positionXMax = max(arrayfun(@(i)max(position{i}(:,1)),1:nStrainlines));
positionYMin = min(arrayfun(@(i)min(position{i}(:,2)),1:nStrainlines));
positionYMax = max(arrayfun(@(i)max(position{i}(:,2)),1:nStrainlines));

output = false;

if positionXMin < domain(1)
    warning([mfilename,':shearlineOutsideDomain'],...
        ['Shearline position outside domain, xMin = ',...
        num2str(positionXMin)]);
    output = true;
end

if positionXMax > domain(3)
    warning([mfilename,':shearlineOutsideDomain'],...
        ['Shearline position outside domain, xMax = ',...
        num2str(positionXMax)]);
    output = true;
end

if positionYMin < domain(2)
    warning([mfilename,':shearlineOutsideDomain'],...
        ['Shearline position outside domain, yMin = ',...
        num2str(positionYMin)]);
    output = true;
end

if positionYMax > domain(4)
    warning([mfilename,':shearlineOutsideDomain'],...
        ['Shearline position outside domain, yMax = ',...
        num2str(positionYMax)]);
    output = true;
end
