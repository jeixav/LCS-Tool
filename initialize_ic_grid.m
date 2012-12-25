%initialize_ic_grid Initialize initial conditions on cartesian grid
function position = initialize_ic_grid(resolution,domain)

xVector = shift_ic(resolution(1),domain(1,:));
yVector = shift_ic(resolution(2),domain(2,:));

[positionX,positionY] = meshgrid(xVector,yVector);
position(:,1) = positionX(:);
position(:,2) = positionY(:);
