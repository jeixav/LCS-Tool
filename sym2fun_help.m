%% sym2fun
% Make a regular MATLAB function from a symbolic function
%
%% Description
% derivative = sym2fun(symDerivative)
%
% Given a symbolic function definition of the flow derivative, create a
% regular function.
%
%% Example
t = sym('t');
x = sym('x');
y = sym('y');
p = struct('epsilon',.1,'a',.1,'omega',pi/5);
forcing = p.epsilon*sin(p.omega*t)*x^2 + (1 - 2*p.epsilon...
   *sin(p.omega*t))*x;
symDerivative(1) = -pi*p.a*sin(pi*forcing)*cos(pi*y);
symDerivative(2) = pi*p.a*cos(pi*forcing).*sin(pi*y)...
   *(2*p.epsilon*sin(p.omega*t)*x + 1 - 2*p.epsilon*sin(p.omega*t));
derivative = sym2fun(symDerivative);
