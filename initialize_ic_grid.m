%initialize_ic_grid Initialize initial conditions on cartesian grid
function position = initialize_ic_grid(resolution,domain)

xVector = linspace(domain(1,1),domain(1,2),resolution(1));
yVector = linspace(domain(2,1),domain(2,2),resolution(2));

[positionX,positionY] = meshgrid(xVector,yVector);
position(:,1) = positionX(:);
position(:,2) = positionY(:);
