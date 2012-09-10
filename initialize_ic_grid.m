function position = initialize_ic_grid(resolution,domain)
%initialize_ic_grid Initialize initial conditions on cartesian grid.

deltaX = (domain(1,2) - domain(1,1))/double(resolution(1));
deltaY = (domain(2,2) - domain(2,1))/double(resolution(2));
[positionX positionY] = meshgrid(...
    domain(1,1) + .5*deltaX:deltaX:domain(1,2) - .5*deltaX,...
    domain(2,1) + .5*deltaY:deltaY:domain(2,2) - .5*deltaY);
position(:,1) = reshape(positionX,numel(positionX),1);
position(:,2) = reshape(positionY,numel(positionY),1);
