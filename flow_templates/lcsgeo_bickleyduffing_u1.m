function lcsgeo_bickleyduffing_u1

% parameters Bickley jet periodic
a = 6371e3;
U = 62.66;
Lx = pi*a;
Ly = 1770e3;
A1 = 0.075;
A2 = 0.4;
A3 = .3;
k1 = 2/a;
k2 = 4/a;
k3 = 6/a;
c2 = .205*U;
c3 = .461*U;
sga1 = .5*k2*(c2-c3);
sga2 = k2*(c2-c3);
T = max(2*pi./abs([sga1 sga2]));

% aperiodic forcing - forced-disipative Duffing oscillator
NT = 5;
w = 5;
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
[tf, xy] = ode23(@(t, xy, dta, bta, alp, gma, oga) ...
   [xy(2); ...
   -dta*xy(2)-bta*xy(1)-alp*xy(1)^3+gma*cos(oga*t)], ...
   0:1e-3:NT*w, [0 0]', ...
   options, ...
   .1, 0, 1, 11, 1);
f = xy(:,2);

% plot
figure(1), clf

% (x,y) for velocity evaluation
x = 2*pi/k1/2;
y = 2.2599*Ly/2;

% appropriate time/amplitude rescaling
t = tf*T/w;
disp(['max(t): ',num2str(max(t))])
F = f/max(f)*1.75;

% strong aperiodic
a2 = A2*F;
a1 = A1*F;
uF = 2*a2*U*sech(y/Ly).^2.*tanh(y/Ly).*cos(k2*x) + ...
   2*a1*U*sech(y/Ly).^2.*tanh(y/Ly).*cos(k1*x);
h1 = plot(t*U/Lx, uF/U, 'LineW', 1.2, 'Color', 'k'); hold on
% weak aperiodic
a2 = A2*F/20;
a1 = A1*F/20;
uF = 2*a2*U*sech(y/Ly).^2.*tanh(y/Ly).*cos(k2*x) + ...
   2*a1*U*sech(y/Ly).^2.*tanh(y/Ly).*cos(k1*x);
h2 = plot(t*U/Lx, uF/U, 'Color', 'k', 'LineS', '--');
% periodic
a2 = A2/10;
a1 = A1/10;
u = 2*a2*U*sech(y/Ly).^2.*tanh(y/Ly).*cos(k2*x-sga2*t) + ...
   2*a1*U*sech(y/Ly).^2.*tanh(y/Ly).*cos(k1*x-sga1*t);
h3 = plot(t*U/Lx, u/U, 'k'); hold off

axis([0 NT*T*U/Lx -inf inf])
box off
xlabel('$tU/L_x$')
ylabel('$u_1^\mathrm{axis}/U$')
legend([h3 h2 h1], {'periodic';'aperiodic (weak)';'aperiodic (strong)';})
