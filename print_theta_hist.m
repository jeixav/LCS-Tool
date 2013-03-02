function print_theta_hist(theta)

[n,xout] = hist(theta,degtorad(5:10:175));

for i = 1:length(n)
    fprintf('%g\t%g\t%g\n',radtodeg(xout(i)),n(i),n(i)/sum(n))
end
