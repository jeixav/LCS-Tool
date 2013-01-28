function travelingWave = traveling_wave
% Produce structures for use with the LCS Tool. Describes a traveling wave
% flow.
%
% Reference: doi:10.1063/1.3272711

flow.parameters = struct('omega',1,...
    'amplitude',1,...
    'waveNumber',1,...
    'speed',.5,...
    'forcingAmplitude',.25);

flow.imposeIncompressibility = false;

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

flow.domain = [-pi pi; 0 pi];
flow.resolution = uint64([2 1]*200);
flow.periodicBc = [true false];

strainline.resolution = uint64([2 1]*5);
strainline.maxLength = 10;
strainline.geodesicDeviationTol = .1;
strainline.lengthTol = 2;
strainline.filteringMethod = 'hausdorff';
strainline.filteringDistanceTol = 1.25;

shearline.resolution = uint64([2 1]);
shearline.finalTime = 2;
shearline.averageGeodesicDeviationNegTol = inf;
shearline.averageGeodesicDeviationPosTol = inf;

noStretchLine.resolution = uint64([2 1]);
noStretchLine.finalTime = 1;

travelingWave.flow = flow;
travelingWave.strainline = strainline;
travelingWave.shearline = shearline;
travelingWave.noStretchLine = noStretchLine;

function symF = symDerivative(parameters)

t = sym('t');
x = sym('x');
y = sym('y');

p = parameters;

forcing = sin(p.omega*t);

symF(1) = p.speed - p.amplitude*sin(p.waveNumber*x)*cos(y) ...
    - p.forcingAmplitude*forcing;
symF(2) = p.amplitude*p.waveNumber*cos(p.waveNumber*x)*sin(y);

function df = symJacDerivative(parameters)

symF = symDerivative(parameters);

x = sym('x');
y = sym('y');

df{1,1} = diff(symF(1),x);
df{1,2} = diff(symF(1),y);
df{2,1} = diff(symF(2),x);
df{2,2} = diff(symF(2),y);
