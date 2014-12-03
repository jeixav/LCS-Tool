% elliptic_lcs Given closed lambda lines, return elliptic LCSs
%
% SYNTAX
% ellipticLcs = elliptic_lcs(closedLambdaLine)
%
% INPUT ARGUMENTS
% closedLambdaLine: closed lambda line positions returned from
% discard_empty_closed_lambda
%
% OUTPUT ARGUMENT
% ellipticLcs: elliptic LCS position for each poincare section over a range
% of lambda values. Cell array with number of elements equal to number of
% Poincare sections
%
% EXAMPLES
% x positions of elliptic LCS for third poincare section:
% [ellipticLcs] = elliptic_lcs(closedLambdaLine)
% ellipticLcs{3}(:,1)
% 
% y positions of elliptic LCS for fourth poincare section:
% ellipticLcs{4}(:,2)
%
% number of elliptic LCSs:
% sum(~cellfun(@isempty,ellipticLcs))
%
% lambda value of closed orbit of second Poincare section:
% lambda(outerPosI(2))

function ellipticLcs = elliptic_lcs(closedLambdaLine)

nLambda = size(closedLambdaLine,1);
nPoincareSection = size(closedLambdaLine,2);

ellipticLcs = cell(1,nPoincareSection);

for iPoincareSection = 1:nPoincareSection
    maxArea = nan;
    for iLambda = 1:nLambda
        outerClosedLambdaLine = closedLambdaLine{iLambda,iPoincareSection}{end};
        iArea = polyarea(outerClosedLambdaLine(:,1),outerClosedLambdaLine(:,2));
        if isnan(maxArea)
            if ~isnan(iArea)
                maxArea = iArea;
                ellipticLcs{iPoincareSection} = outerClosedLambdaLine;
            end
        else
            if ~isnan(iArea)
                if iArea > maxArea
                    maxArea = iArea;
                    ellipticLcs{iPoincareSection} = outerClosedLambdaLine;
                end
            end
        end
    end
end

% Discard ellipticLcs elements for Poincare sections without closed orbits
ellipticLcs = ellipticLcs(~cellfun(@isempty,ellipticLcs));
