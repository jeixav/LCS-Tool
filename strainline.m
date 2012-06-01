function strainline(flow,position,eigenvalues,eigenvectors)

% FIXME How is timestep chosen
timestep = .01;
tspan = flow.initial_time:timestep:flow.final_time;

n_x = strainlines.resolution(1);
n_y = strainlines.resolution(2);
[px py] = meshgrid(linspace(flow.domain(1,1),flow.domain(1,2),n_x),...
    linspace(flow.domain(2,1),flow.domain(2,2),n_y));
px = px(2:end-1,2:end-1);
py = py(2:end-1,2:end-1);
px = reshape(px,numel(px),1);
py = reshape(py,numel(py),1);
    
[x_fw y_fw] = xi1_tracing([px py],tspan,position,flow.resolution,...
    eigenvalues,eigenvectors);
[x_bw y_bw] = xi1_tracing([px py],-tspan,position,flow.resolution,...
    eigenvalues,eigenvectors);
% [x_bw y_bw] = xi1_tracing(px,py,-tspan,xi,yi,l1,l2,v1);

% Discard first line of backward time since it repeats forward time
if any(x_fw(1,:) ~= x_bw(1,:)) || any(y_fw(1,:) ~= y_bw(1,:))
    warning('Main:FirstLinesUnequal', ...
        'Initial positions of forward and backward times unequal')
else
    x_st = ([flipud(x_fw); x_bw(2:end,:)]);
    y_st = ([flipud(y_fw); y_bw(2:end,:)]);
end

figure
a1 = axes;
set(a1,'nextplot','add')
plot_initial_conditions = true;
if plot_initial_conditions
    p1 = plot(px,py);
    set(p1,'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'Marker', 'o', ...
        'LineStyle', 'none')
end

plot(x_st,y_st,'k');

end
