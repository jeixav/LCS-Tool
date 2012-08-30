function doubleGyre = double_gyre
% Describes a double gyre flow.
%
% References: doi:10.1016/j.physd.2005.10.007, doi:10.5194/npg-7-59-2000,
% doi:10.5194/npg-4-223-1997

flow.parameters = struct('epsilon',.1,...
    'a',.1,...
    'omega',pi/5);

flow.isCompressible = false;

symF = symDerivative(flow.parameters);
dyScalar1 = matlabFunction(symF(1),'vars',{'t','x','y'});
dyScalar2 = matlabFunction(symF(2),'vars',{'t','x','y'});

flow.derivative = @(t,y)[dyScalar1(t,y(1),y(2)); dyScalar2(t,y(1),y(2))];

symJacDy = symJacDerivative(flow.parameters);

jacDyScalar11 = matlabFunction(symJacDy{1,1},'vars',{'t','x','y'});
jacDyScalar12 = matlabFunction(symJacDy{1,2},'vars',{'t','x','y'});
jacDyScalar21 = matlabFunction(symJacDy{2,1},'vars',{'t','x','y'});
jacDyScalar22 = matlabFunction(symJacDy{2,2},'vars',{'t','x','y'});

flow.dDerivative = @(t,y)[jacDyScalar11(t,y(1),y(2)) ...
    jacDyScalar12(t,y(1),y(2)); jacDyScalar21(t,y(1),y(2)) ...
    jacDyScalar22(t,y(1),y(2))];

flow.timespan = [0 20];

flow.domain = [0 2; 0 1];
flow.resolution = uint64([2 1]*20);

strainline.resolution = uint64([2 1]*5);
strainline.finalTime = 10;
strainline.geodesicDeviationTol = inf;
strainline.lengthTol = 0;
strainline.filteringMethod = 'hausdorff';
strainline.filteringDistanceTol = 0;
strainline.odeSolverOptions = odeset('relTol',1e-4);

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

function symF = symDerivative(parameters)

t = sym('t');
x = sym('x');
y = sym('y');

p = parameters;

forcing = p.epsilon*sin(p.omega*t)*x^2 + (1 - 2*p.epsilon...
    *sin(p.omega*t))*x;

symF(1) = -pi*p.a*sin(pi*forcing)*cos(pi*y);
symF(2) = pi*p.a*cos(pi*forcing).*sin(pi*y)...
    *(2*p.epsilon*sin(p.omega*t)*x + 1 - 2*p.epsilon*sin(p.omega*t));

function df = symJacDerivative(parameters)

symF = symDerivative(parameters);

x = sym('x');
y = sym('y');

df{1,1} = diff(symF(1),x);
df{1,2} = diff(symF(1),y);
df{2,1} = diff(symF(2),x);
df{2,2} = diff(symF(2),y);
