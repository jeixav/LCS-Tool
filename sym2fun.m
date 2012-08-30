% sym2fun Make a regular MATLAB function from a symbolic function
%
% DESCRIPTION
% ...
%
% EXAMPLE
% ...
%

function derivative = sym2fun(symDerivative)

dyScalar1 = matlabFunction(symDerivative(1),'vars',{'t','x','y'});
dyScalar2 = matlabFunction(symDerivative(2),'vars',{'t','x','y'});

derivative = @(t,y)[dyScalar1(t,y(1),y(2)); dyScalar2(t,y(1),y(2))];
