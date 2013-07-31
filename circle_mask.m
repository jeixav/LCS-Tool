function circle = circle_mask(radius)

[X,Y] = meshgrid(0:radius,radius:-1:0);
distance = sqrt(double(X).^2 + double(Y).^2);
quarterCircle = distance <= radius;
circle = [[fliplr(quarterCircle),quarterCircle(:,2:end)];[rot90(quarterCircle(1:end-1,:),2),flipud(quarterCircle(1:end-1,2:end))]];
