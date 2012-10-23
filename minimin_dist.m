% Calculates the minimin distance between P and Q

function d = minimin_dist(p,q)

sP = size(p);
sQ = size(q);

if ~(sP(2)==sQ(2))
    error('Inputs P and Q must have the same number of columns')
end

% obtain all possible point comparisons
iP = repmat(1:sP(1),[1,sQ(1)])';
iQ = repmat(1:sQ(1),[sP(1),1]);
combos = [iP,iQ(:)];

% get distances for each point combination
cP = p(combos(:,1),:);
cQ = q(combos(:,2),:);
dists = sqrt(sum((cP - cQ).^2,2));

% Now create a matrix of distances where D(n,m) is the distance of the nth
% point in P from the mth point in Q. The maximum distance from any point
% in Q from P will be max(D,[],1) and the maximum distance from any point
% in P from Q will be max(D,[],1);
d = reshape(dists,sP(1),[]);

d = min(d(:));
