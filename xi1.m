function dy = xi1(~,y,v1,xi,yi)

global XI1_1B XI1_2B

N = round(length(y)/2);

% freeze the particles at the boundary
y(1:N,1) = min(y(1:N,1),max(max(xi)));
y(1:N,1) = max(y(1:N,1),min(min(xi)));
y(N+1:2*N,1) = min(y(N+1:2*N,1),max(max(yi)));
y(N+1:2*N,1) = max(y(N+1:2*N,1),min(min(yi)));

dy = zeros(2*N,1);

% xi1_1 = interp2(xi,yi,v1(:,:,1),y(1:N,1),y(N+1:2*N,1),interp_method);
% xi1_2 = interp2(xi,yi,v1(:,:,2),y(1:N,1),y(N+1:2*N,1),interp_method);
[xi1_1 xi1_2] = smooth_xi1(y(1:N,1),y(N+1:2*N,1),xi,yi,...
    cat(3,v1(:,:,1),v1(:,:,2)));

dy(1:N,1) = sign(xi1_1.*XI1_1B + xi1_2.*XI1_2B).*xi1_1;
dy(N+1:2*N,1) = sign(xi1_1.*XI1_1B + xi1_2.*XI1_2B).*xi1_2;

XI1_1B = sign(xi1_1.*XI1_1B + xi1_2.*XI1_2B).*xi1_1;
XI1_2B = sign(xi1_1.*XI1_1B + xi1_2.*XI1_2B).*xi1_2;
