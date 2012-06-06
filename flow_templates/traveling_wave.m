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

flow.isCompressible = false;

flow.derivative = @derivative;

flow.timespan = [0 20];

flow.domain = [-pi pi; 0 pi];
flow.resolution = uint64([2 1]*200);
flow.periodicBc = [true false];

flow.auxiliaryGridRelativeDelta = 5e-3;

strainline.resolution = uint64([2 1]*5);
% strainline.timestep = .025;
strainline.finalTime = 10;
strainline.geodesicDeviationTol = .1;
strainline.lengthTol = 2;
strainline.filteringMethod = 'hausdorff';
strainline.filteringDistanceTol = 1.25;

shearline.resolution = uint64([2 1]);
shearline.finalTime = 10;
shearline.averageGeodesicDeviationNegTol = inf;
shearline.averageGeodesicDeviationPosTol = inf;

travelingWave.flow = flow;
travelingWave.strainline = strainline;
travelingWave.shearline = shearline;

function dy = derivative(t,y,parameters)

N = length(y)/2;
x1 = y(1:N,1);
x2 = y(N+1:2*N,1);

omega = parameters.omega;
amplitude = parameters.amplitude;
waveNumber = parameters.waveNumber;
speed = parameters.speed;
forcingAmplitude = parameters.forcingAmplitude;

forcing = sin(omega*t);

dy = nan(2*N,1);
dy(1:N,1) = speed - amplitude*sin(waveNumber*x1).*cos(x2) ...
    - forcingAmplitude*forcing;
dy(N+1:2*N,1) = amplitude*waveNumber*cos(waveNumber*x1).*sin(x2);

% x1 = y(1);
% x2 = y(2);
% 
% omega = parameters.omega;
% amplitude = parameters.amplitude;
% waveNumber = parameters.waveNumber;
% speed = parameters.speed;
% forcingAmplitude = parameters.forcingAmplitude;
% 
% forcing = sin(omega*t);
% 
% dy(1) = speed - amplitude*sin(waveNumber*x1)*cos(x2) ...
%     - forcingAmplitude*forcing;
% dy(2) = amplitude*waveNumber*cos(waveNumber*x1)*sin(x2);
