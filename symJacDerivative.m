function df = symJacDerivative(symDerivative)

x = sym('x');
y = sym('y');

df{1,1} = diff(symDerivative(1),x);
df{1,2} = diff(symDerivative(1),y);
df{2,1} = diff(symDerivative(2),x);
df{2,2} = diff(symDerivative(2),y);
