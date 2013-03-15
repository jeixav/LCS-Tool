% Cauchy-Green strain tensor statistics

function o = cgStrain_stats(cgStrain,cgStrainEigenvector,cgStrainEigenvalue,verbose)

%% Negative eigenvalues
if verbose
    fprintf('Number of negative eigenvalues: %u.\n',numel(find(cgStrainEigenvalue(:) < 0)))
end

%% Eigenvalue and eigenvector error
n = size(cgStrain,3);

    function CGEigErrorArrayfun = cg_eig_error_arrayfun(idx,cgStrain,cgStrainEigenvector,cgStrainEigenvalue)
        CGEigErrorArrayfun = eig_error(cgStrain(:,:,idx),reshape(cgStrainEigenvector(idx,:),2,2),diag(cgStrainEigenvalue(idx,:)));
    end

EigError = arrayfun(@(idx)cg_eig_error_arrayfun(idx,cgStrain,cgStrainEigenvector,cgStrainEigenvalue),1:n,'UniformOutput',false);
EigError = cell2mat(transpose(EigError));

if verbose
    fprintf('max(cgStrain*cgStrainEigenvector - cgStrainEigenvalue*cgStrainEigenvector):\n')
    fprintf('\t1: %g\n',max(EigError(:,1)))
    fprintf('\t2: %g\n',max(EigError(:,2)))
end

%% Product of eigenvalues
prodCgStrainD = prod(cgStrainEigenvalue,2);

o(1) = 1 - max(cgStrainEigenvalue(:,1));
detCgStrain = arrayfun(@(idx)det(cgStrain(:,:,idx)),1:size(cgStrain,3));
o(2) = mean(abs(detCgStrain-1));

if verbose
    fprintf('lambda_1*lambda2:\n')
    fprintf('\tmin = %g\n',min(prodCgStrainD))
    fprintf('\tmax = %g\n',max(prodCgStrainD))
    fprintf('\tmean = %g\n',mean(prodCgStrainD))
    fprintf('\tmedian = %g\n',median(prodCgStrainD))
    fprintf('1 - max(lambda_1) = %g\n',o(1))
    fprintf('mean(abs(detCgStrain-1)) = %g\n',o(2))
end

end
