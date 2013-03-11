% error Error measurement of eigenvalues and eigenvectors
%
% SYNTAX
% e = eig_error(a,v,d)
%
% DESCRIPTION
% The error is defined by:
% e(i) = norm((a - d(i,i)*eye(size(a)))*v(:,i))
%
% EXAMPLE
% a = rand(4);
% [v,d] = eig(a);
% e = eig_error(a,v,d);

function EigError = eig_error(a,v,d)

    function EigErrorArrayfun = eig_error_arrayfun(idx,a,v,d)
        EigErrorArrayfun = norm((a - d(idx,idx)*eye(size(a)))*v(:,idx));
    end        

EigError = arrayfun(@(idx)eig_error_arrayfun(idx,a,v,d),1:size(a,1));

end
