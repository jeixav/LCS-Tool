% shift_ic Place all positions inside domain
function vector = shift_ic(resolution,domain)

delta = (domain(2) - domain(1))/double(resolution);

vector = linspace(domain(1),domain(2),resolution+1);
vector = vector + .5*delta;
vector = vector(1:end-1);