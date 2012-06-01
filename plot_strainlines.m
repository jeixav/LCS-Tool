function plot_strainlines(axes,position)

plot(axes,position(:,:,1),position(:,:,1),'k','Tag','strainline');

plot_initial_conditions = false;
if plot_initial_conditions
    p1 = plot(axes,px,py);
    set(p1,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','o',...
        'LineStyle','none')
end

end
