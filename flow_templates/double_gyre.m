function doubleGyre = double_gyre
% Describes a double gyre flow.
%
% References: doi:10.1016/j.physd.2005.10.007, doi:10.5194/npg-7-59-2000,
% doi:10.5194/npg-4-223-1997

flow.parameters = struct('epsilon',.1,...
    'a',.1,...
    'omega',2*pi/10);

flow.isCompressible = false;

flow.timespan = [0 20];
flow.derivative = @derivative;

flow.domain = [0 2; 0 1];
flow.resolution = uint64([2 1]*200);

flow.auxiliaryGridRelativeDelta = 5e-3;

strainline.resolution = uint64([2 1]*5);
strainline.finalTime = 10;
strainline.geodesicDeviationTol = .1;
strainline.lengthTol = 2;
strainline.filteringMethod = 'hausdorff';
strainline.filteringDistanceTol = 1.25;

shearline.resolution = uint64([2 1]);
shearline.finalTime = 1;
shearline.averageGeodesicDeviationNegTol = inf;
shearline.averageGeodesicDeviationPosTol = inf;

noStretchLine.resolution = uint64([2 1]);
noStretchLine.finalTime = 1;

doubleGyre.flow = flow;
doubleGyre.strainline = strainline;
doubleGyre.shearline = shearline;
doubleGyre.noStretchLine = noStretchLine;

function dy = derivative(t,y,parameters)

N = length(y)/2;
x1 = y(1:N,1);
x2 = y(N+1:2*N,1);

epsilon = parameters.epsilon;
a = parameters.a;
omega = parameters.omega;

forcing = epsilon*sin(omega*t)*x1.^2 + (1 - 2*epsilon*sin(omega*t))*x1;

dy = nan(2*N,1);
dy(1:N,1) = -pi*a*sin(pi*forcing).*cos(pi*x2);
dy(N+1:2*N,1) = pi*a*cos(pi*forcing).*sin(pi*x2)...
    .*(2*epsilon*sin(omega*t)*x1 + 1 - 2*epsilon*sin(omega*t));
