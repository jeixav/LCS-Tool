function supermin_demo(varargin)
%SUPERMIN_DEMO Demonstrate superminimization graphically.

switch nargin
    case 0
        superminDistance = 1;
    case 1
        superminDistance = varargin{1};
    otherwise
        error('Incorrect number of input arguments')
end

showSuperminDataPlot = true;

[superminData,a1] = supermin_data(showSuperminDataPlot);

superminIndex = superminimize(superminData.position,...
    superminData.value,superminDistance);

plot(a1,superminData.position(superminIndex),...
    superminData.value(superminIndex),...
    'linestyle','none','marker','o','markerfacecolor','r')

end

