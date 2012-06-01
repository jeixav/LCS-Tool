function length = curve_length(position)

x_d = diff(position(:,1));
y_d = diff(position(:,2));

length = sum(sqrt(x_d.^2 + y_d.^2));

end
