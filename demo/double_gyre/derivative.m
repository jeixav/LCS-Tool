% derivative Double gyre velocity field
%
% SYNTAX
% derivative_ = derivative(t,position,useEoV,epsilon,amplitude,omega)
%
% INPUT ARGUMENTS
% t: time
% position: [x1;y1;x2;y2;...;xn;yn]
% useEov: logical that controls use of the equation of variation
% epsilon,amplitude,omega: double gyre parameters
%
% REFERENCE
% DOI:10.1016/j.physd.2005.10.007

function derivative_ = derivative(t,x,useEoV,epsilon,amplitude,omega)

validateattributes(t,{'double'},{'scalar'})
validateattributes(x,{'double'},{'column'})
% Cannot use validateattributes to check x is a column vector with an even
% number of elements
if mod(numel(x),2)
    error([mfilename,':numelXNotEven'],'numel(x) = %d is not even',numel(x))
end
validateattributes(x,{'double'},{'column'})
validateattributes(useEoV,{'logical'},{'scalar'})
validateattributes(epsilon,{'double'},{'scalar'})
validateattributes(amplitude,{'double'},{'scalar'})
validateattributes(omega,{'double'},{'scalar'})

if useEoV
    idx1 = 1:6:size(x,1)-5;
    idx2 = 2:6:size(x,1)-4;
else
    idx1 = 1:2:numel(x)-1;
    idx2 = 2:2:numel(x);
end

a = epsilon*sin(omega*t);
b = 1 - 2*epsilon*sin(omega*t);
forcing = a*x(idx1).^2 + b*x(idx1);

derivative_ = nan(size(x));

derivative_(idx1) = -pi*amplitude*sin(pi*forcing).*cos(pi*x(idx2));
derivative_(idx2) = pi*amplitude*cos(pi*forcing).*sin(pi*x(idx2)).*(2*a*x(idx1) + b);

if useEoV
    % Define terms of the equation of variation
    idx3 = 3:6:size(x,1)-3;
    idx4 = 4:6:size(x,1)-2;
    idx5 = 5:6:size(x,1)-1;
    idx6 = 6:6:size(x,1);
    
    dux = -pi^2*amplitude*cos(pi*forcing).*cos(pi*x(idx2)).*(2*a*x(idx1) + b);
    duy = pi^2*amplitude*sin(pi*forcing).*sin(pi*x(idx2));
    dvx = -pi^2*amplitude*sin(pi*forcing).*sin(pi*x(idx2)).*(2*a*x(idx1) + b) + 2*a*pi*amplitude*cos(pi*forcing).*sin(pi*x(idx2));
    dvy = pi^2*amplitude*cos(pi*forcing).*cos(pi*x(idx2)).*(2*a*x(idx1) + b);
    
    % Perform matrix multiplication manually
    derivative_(idx3) = dux.*x(idx3) + duy.*x(idx5);
    derivative_(idx4) = dux.*x(idx4) + duy.*x(idx6);
    derivative_(idx5) = dvx.*x(idx3) + dvy.*x(idx5);
    derivative_(idx6) = dvx.*x(idx4) + dvy.*x(idx6);
end
