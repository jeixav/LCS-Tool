function row_array = from_coupled(coupled_array,auxiliary_array_flag)
% FROM_COUPLED Reshape positions array from coupled ODE form
%    ROW_ARRAY = FROM_COUPLED(COUPLED_ARRAY,AUXILIARY_ARRAY_FLAG)
%
%    Set AUXILIARY_ARRAY_FLAG = TRUE if COUPLED_ARRAY represents an array
%    with auxiliary grid points.

dimensions = 2;
if auxiliary_array_flag
    num_aux_coords = 4;
else
    num_aux_coords = 1;
end

num = size(coupled_array,1)/num_aux_coords/dimensions;

x = transpose(reshape(coupled_array(1:end/2),num_aux_coords,num));
y = transpose(reshape(coupled_array(end/2+1:end),num_aux_coords,num));
clear('coupled_array')

row_array = nan(num,num_aux_coords*dimensions);
row_array(:,1:2:end-1) = x;
row_array(:,2:2:end) = y;

end
