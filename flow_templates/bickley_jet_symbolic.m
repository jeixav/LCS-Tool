%bickley_jet_symbolic Helper function to create Bickley jet derivatives

function bickley_jet_symbolic

syms xVar y c_3 U L_y epsilon_1 epsilon_2 epsilon_3 k_1 k_2 k_3 f_1 f_2

psi_0 = c_3*y - U*L_y*tanh(y/L_y) + epsilon_3*U*L_y*sech(y/L_y)*cos(k_3*xVar);

psi_1 = U*L_y*sech(y/L_y)^2*real(epsilon_1*f_1*exp(1i*k_1*xVar) + epsilon_2*f_2*exp(1i*k_2*xVar));

psi = psi_0 + psi_1;

u = simple(diff(-psi,y));
v = diff(psi,xVar);

u_x = diff(u,xVar);
u_y = diff(u,y);
v_x = diff(v,xVar);
v_y = diff(v,y);

disp('u:')
uStr = all_strrep(char(u));
disp(uStr)

disp('v:')
vStr = all_strrep(char(v));
disp(vStr)

disp('u_x:')
u_xStr = all_strrep(char(u_x));
disp(u_xStr)

disp('u_y:')
u_yStr = all_strrep(char(u_y));
disp(u_yStr)

disp('v_x:')
v_xStr = all_strrep(char(v_x));
disp(v_xStr)

disp('v_y:')
v_yStr = all_strrep(char(v_y));
disp(v_yStr)

function modifiedStr = all_strrep(origStr)

modifiedStr = strrep(origStr,'U','p.u');
modifiedStr = strrep(modifiedStr,'L_y','p.lengthY');
modifiedStr = strrep(modifiedStr,'c_3','p.c3');
modifiedStr = strrep(modifiedStr,'epsilon_1','p.epsilon1');
modifiedStr = strrep(modifiedStr,'epsilon_2','p.epsilon2');
modifiedStr = strrep(modifiedStr,'epsilon_3','p.epsilon3');
modifiedStr = strrep(modifiedStr,'k_1','k(1)');
modifiedStr = strrep(modifiedStr,'k_2','k(2)');
modifiedStr = strrep(modifiedStr,'k_3','k(3)');
modifiedStr = strrep(modifiedStr,'f_1','f1(t)');
modifiedStr = strrep(modifiedStr,'f_2','f2(t)');
modifiedStr = strrep(modifiedStr,'xVar','x(idx1)');
modifiedStr = strrep(modifiedStr,'y','x(idx2)');

modifiedStr = strrep(modifiedStr,'*','.*');
modifiedStr = strrep(modifiedStr,'^','.^');
modifiedStr = strrep(modifiedStr,'/','./');
