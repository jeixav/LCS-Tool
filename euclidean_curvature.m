% euclidean_curvature Euclidean curvature of a trajectory of the vector2
% vector field
%
% SYNTAX
% kappa = euclidean_curvature(vector1,vector2,spacing)
%
% DESCRIPTION
% Equation (30) from doi:10.1016/j.physd.2012.06.012. Orientation
% discontinuities in the eigenvector field are taken into account.

function kappa = euclidean_curvature(vector1,vector2,spacing)

rotation_check(vector1,vector2)

resolution(1) = size(vector1,3);
resolution(2) = size(vector1,2);

gradCgEigenvector1 = nan(2,2,resolution(2),resolution(1));

% x-direction gradient
ref = vector1(:,:,2:end-1);

refP = vector1(:,:,3:end);
s = sign(squeeze(dot(ref,refP)));
s = shiftdim(repmat(s,[1,1,2]),2);
refP = s.*refP;

refM = vector1(:,:,1:end-2);
s = sign(squeeze(dot(ref,refM)));
s = shiftdim(repmat(s,[1,1,2]),2);
refM = s.*refM;

gradCgEigenvector1(:,1,:,2:end-1) = .5*(refP - refM)/spacing(1);

% y-direction gradient
ref = vector1(:,2:end-1,:);

refP = vector1(:,3:end,:);
s = sign(squeeze(dot(ref,refP)));
s = shiftdim(repmat(s,[1,1,2]),2);
refP = s.*refP;

refM = vector1(:,1:end-2,:);
s = sign(squeeze(dot(ref,refM)));
s = shiftdim(repmat(s,[1,1,2]),2);
refM = s.*refM;

gradCgEigenvector1(:,2,2:end-1,:) = .5*(refP - refM)/spacing(2);

% One-sided finite-difference on boundaries
for m = [1,resolution(1)]
    ref = vector1(:,:,m);

    if m == 1
        refP = vector1(:,:,m+1);
        s = sign(squeeze(dot(ref,refP)));
    else
        refM = vector1(:,:,m-1);
        s = sign(squeeze(dot(ref,refM)));
    end
    s = repmat(s,[2,1]);

    if m == 1
        refP = s.*refP;
        gradCgEigenvector1(:,1,:,m) = (refP-ref)/spacing(1);
    else
        refM = s.*refM;
        gradCgEigenvector1(:,1,:,m) = (ref-refM)/spacing(1);
    end
end

for n = [1,resolution(2)]
    ref = squeeze(vector1(:,n,:));
    
    if n == 1
        refP = squeeze(vector1(:,n+1,:));
        s = sign(squeeze(dot(ref,refP)));
    else
        refM = squeeze(vector1(:,n-1,:));
        s = sign(squeeze(dot(ref,refM)));
    end
    s = repmat(s,[2,1]);

    if n == 1
        refP = s.*refP;
        gradCgEigenvector1(:,2,n,:) = (refP-ref)/spacing(2);
    else
        refM = s.*refM;
        gradCgEigenvector1(:,2,n,:) = (ref-refM)/spacing(2);
    end
end

% One-sided differences at corners
% Corner(1,1)
m = 1;
n = 1;

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m+1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (b(1)-a(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (b(2)-a(2))/spacing(1);

b = squeeze(vector1(:,n+1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (b(1)-a(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (b(2)-a(2))/spacing(2);

% Corner(1,res(2))
n = resolution(2);

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m+1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (b(1)-a(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (b(2)-a(2))/spacing(1);

b = squeeze(vector1(:,n-1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (a(1)-b(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (a(2)-b(2))/spacing(2);

% Corner(res(1),res(2))
m = resolution(1);

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m-1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (a(1)-b(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (a(2)-b(2))/spacing(1);

b = squeeze(vector1(:,n-1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (a(1)-b(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (a(2)-b(2))/spacing(2);

% Corner(res(1),1)
n = 1;

a = squeeze(vector1(:,n,m));

b = squeeze(vector1(:,n,m-1));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,1,n,m) = (a(1)-b(1))/spacing(1);
gradCgEigenvector1(2,1,n,m) = (a(2)-b(2))/spacing(1);

b = squeeze(vector1(:,n+1,m));
if sign(dot(a,b)) < 0
    b = -b;
end

gradCgEigenvector1(1,2,n,m) = (b(1)-a(1))/spacing(2);
gradCgEigenvector1(2,2,n,m) = (b(2)-a(2))/spacing(2);

% Equation 30 from doi:10.1016/j.physd.2012.06.012
kappa = (squeeze(gradCgEigenvector1(1,1,:,:)).*squeeze(vector1(1,:,:)) + squeeze(gradCgEigenvector1(1,2,:,:)).*squeeze(vector1(2,:,:))).*squeeze(vector2(1,:,:)) + (squeeze(gradCgEigenvector1(2,1,:,:)).*squeeze(vector1(1,:,:)) + squeeze(gradCgEigenvector1(2,2,:,:)).*squeeze(vector1(2,:,:))).*squeeze(vector2(2,:,:));

% Ensure when vector1 has an orientation discontinuity, so does vector2.
% Check this by verifying vector2 is consistently rotated in the same
% direction (by ±90°) compared to vector1.
%
% The custom method used to calculate eigenvalues should guarantee this
% always. MATLAB's eig function seems to satisfy this condition as well,
% but it is not documented, so it is best to verify.
function rotation_check(vector1,vector2)

% Set rotation based on that of first grid point vectors
rotationMatrix = [0,1;-1,0];
if any(rotationMatrix*vector2(:,1) ~= vector1(:,1))
    rotationMatrix = [0,-1;1,0];
    if any(rotationMatrix*vector2(:,1) ~= vector1(:,1))
        error(['Unable to determine 90° rotation direction from vector1(:,1) = ',num2str(transpose(vector1(:,1))),' and vector1(:,2) = ',num2str(transpose(vector2(:,1))),'.'])
    end
end

nanIdx = isnan(vector1(1,:));

if any(any(rotationMatrix*vector2(:,~nanIdx) ~= vector1(:,~nanIdx)))
    error('vector1 not rotated by 90° with respect to vector2')
end
