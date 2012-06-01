function coupled_array = to_coupled(row_array)
% TO_COUPLED Reshape positions array to coupled ODE form
%    COUPLED_ARRAY = TO_COUPLED(ROW_ARRAY)

aux_pos_x =  transpose(row_array(:,1:2:end-1));
aux_pos_y =  transpose(row_array(:,2:2:end));
clear('row_array')

coupled_array = [aux_pos_x(:); aux_pos_y(:)];
clear('aux_pos_x','aux_pos_y')

end
