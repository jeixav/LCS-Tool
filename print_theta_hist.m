% print_theta_hist Print a list of angular change values
%
% SYNTAX
% print_theta_hist(theta)
%
% EXAMPLE
% bickleyJet = bickley_jet_coupled(3);
% fullDomain = bickleyJet.flow.domain;
% domain = [2.6,2.7;-.04,.01];
% bickleyJet.flow = set_flow_domain(domain,bickleyJet.flow);
% bickleyJet.flow = set_flow_resolution(64,bickleyJet.flow);
% bickleyJet.strainline = set_strainline_resolution([5 4],bickleyJet.strainline);
% showPlot.quiver = true;
% showPlot.strainline = true;
% [bickleyJet,hAxes] = strain_lcs_script(bickleyJet,showPlot);
% set(findobj(hAxes,'tag','strainline'),'color','k')
% xi1 = bickleyJet.flow.cgEigenvector(:,[1,2]);
% xi1 = shiftdim(reshape(xi1,[fliplr(bickleyJet.flow.resolution),2]),2);
% [thetaX,thetaY] = angle_change(xi1);
% print_theta_hist([thetaX(:);thetaY(:)])

function print_theta_hist(theta)

[n,xout] = hist(theta,degtorad(5:10:175));

for i = 1:length(n)
    fprintf('%g\t%g\t%g\n',radtodeg(xout(i)),n(i),n(i)/sum(n))
end
