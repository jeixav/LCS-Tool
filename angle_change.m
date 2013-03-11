% angle_change Angular change in vector field at adjacent grid points.
%
% SYNTAX
% theta = angle_change(vector)
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
% xi1 = reshape(xi1,[fliplr(bickleyJet.flow.resolution),2]);
% theta = angle_change(xi1);
% print_theta_hist(theta)

function theta = angle_change(vector)

% Compare x to x + DeltaX
a = vector(:,1:end-1,:);
b = vector(:,2:end,:);

normA = norm_array(a);
normB = norm_array(b);

% Formula used is: dot(a,b) = norm(a)*norm(b)*cos(theta)
thetaX = acos(dot(a,b,3)./normA./normB);

% Compare y to y + DeltaY
a = vector(1:end-1,:,:);
b = vector(2:end,:,:);

normA = norm_array(a);
normB = norm_array(b);

thetaY = acos(dot(a,b,3)./normA./normB);

theta = [thetaX(:);thetaY(:)];

% Norm of m-by-n-by-2 array.
function normArray = norm_array(array)

normArray = nan(size(array,1),size(array,2));

% FIXME Using arrayfun would probably be faster than for loops
for idx1 = 1:size(array,1)
    for idx2 = 1:size(array,2)
        normArray(idx1,idx2) = norm(squeeze(array(idx1,idx2,:)));
    end
end
