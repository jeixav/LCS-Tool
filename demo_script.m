doubleGyre = double_gyre;
doubleGyre = set_flow_resolution(uint64([2 1]*500),doubleGyre);
doubleGyre.shearline = set_shearline_ode_solver_options(odeset('relTol',...
    1e-8,'absTol',1e-8),doubleGyre.shearline);
doubleGyre.shearline = set_shearline_initial_position([1+1/3 .5; 1.4 .3;...
    1.5 .4; 1.5 .225; 1.5 .65],doubleGyre.shearline);
doubleGyre = shear_lcs_script(doubleGyre);
set(findobj(gca,'tag','shearlinePosFiltered'),'visible','off')
plot(doubleGyre.shearline.initialPosition(:,1),...
    doubleGyre.shearline.initialPosition(:,2),'ko','tag','initialPosition')