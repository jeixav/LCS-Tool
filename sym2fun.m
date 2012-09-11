% sym2fun Make a regular MATLAB function from a symbolic function
%
% DESCRIPTION
% derivative = sym2fun(symDerivative)
%
% Given a symbolic function definition of the flow derivative, create a
% regular function.
%
% EXAMPLE
% t = sym('t');
% x = sym('x');
% y = sym('y');
% p = struct('epsilon',.1,'a',.1,'omega',pi/5);
% forcing = p.epsilon*sin(p.omega*t)*x^2 + (1 - 2*p.epsilon...
%    *sin(p.omega*t))*x;
% symDerivative(1) = -pi*p.a*sin(pi*forcing)*cos(pi*y);
% symDerivative(2) = pi*p.a*cos(pi*forcing).*sin(pi*y)...
%    *(2*p.epsilon*sin(p.omega*t)*x + 1 - 2*p.epsilon*sin(p.omega*t));
% derivative = sym2fun(symDerivative);

function derivative = sym2fun(symDerivative)

dyScalar1 = matlabFunction(symDerivative(1),'vars',{'t','x','y'});
dyScalar2 = matlabFunction(symDerivative(2),'vars',{'t','x','y'});

derivative = @(t,y)[dyScalar1(t,y(1),y(2)); dyScalar2(t,y(1),y(2))];
