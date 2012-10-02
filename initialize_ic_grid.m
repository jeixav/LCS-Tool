%initialize_ic_grid Initialize initial conditions on cartesian grid
function position = initialize_ic_grid(resolution,domain)

xVector = shift_ic(resolution(1),domain(1,:));
yVector = shift_ic(resolution(2),domain(2,:));

[positionX positionY] = meshgrid(xVector,yVector);
position(:,1) = positionX(:);
position(:,2) = positionY(:);

% shift_ic Place all positions inside domain
function vector = shift_ic(resolution,domain)

delta = (domain(2) - domain(1))/double(resolution);

vector = linspace(domain(1),domain(2),resolution+1);
vector = vector + .5*delta;
vector = vector(1:end-1);