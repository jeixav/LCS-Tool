function  dy = flowdata_derivative(t,y, dummy, VX_interpolant, VY_interpolant)
N = length(y);
dy = nan(N,1);

%% positions are integrated at once (no parfor loop)
% lon, x
dy(1:2:N-1) = VX_interpolant(t*ones(N/2,1),y(2:2:N),y(1:2:N-1));
% lat, y
dy(2:2:N) = VY_interpolant(t*ones(N/2,1),y(2:2:N),y(1:2:N-1));

end