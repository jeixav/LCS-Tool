% plot_elliptic_lcs Plot elliptic LCSs
%
% SYNTAX
% h = plot_elliptic_lcs(hAxes,ellipticLcs)
%
% INPUT ARGUMENTS
% ellipticLcs: elliptic LCS positions returned from elliptic_lcs

function h = plot_elliptic_lcs(hAxes,ellipticLcs)

nPoincareSection = numel(ellipticLcs);

h = gobjects(1,nPoincareSection);

for iPs = 1:nPoincareSection
    h(iPs) = plot(hAxes,ellipticLcs{iPs}(:,1),ellipticLcs{iPs}(:,2));
end

h = h(isgraphics(h));
