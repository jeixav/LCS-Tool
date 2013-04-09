% cgStrain_stats Cauchy-Green strain tensor statistics
%
% SYNTAX
% cgStrain_stats(cgStrain,cgStrainEigenvector,cgStrainEigenvalue)

function cgStrain_stats(cgStrain,cgStrainEigenvector,cgStrainEigenvalue)

%% Negative eigenvalues
fprintf('Number of negative eigenvalues: %u.\n',numel(find(cgStrainEigenvalue(:) < 0)))

nanIdx = isnan(cgStrainEigenvalue(:,1));
fprintf('Number of NaN values: %u.\n',sum(nanIdx))

%% Eigenvalue and eigenvector error
n = size(cgStrain,3);

    function CGEigErrorArrayfun = cg_eig_error_arrayfun(idx,cgStrain,cgStrainEigenvector,cgStrainEigenvalue)
        CGEigErrorArrayfun = eig_error(cgStrain(:,:,idx),reshape(cgStrainEigenvector(idx,:),2,2),diag(cgStrainEigenvalue(idx,:)));
    end

EigError = arrayfun(@(idx)cg_eig_error_arrayfun(idx,cgStrain,cgStrainEigenvector,cgStrainEigenvalue),1:n,'UniformOutput',false);
EigError = cell2mat(transpose(EigError));

fprintf('cgStrain*cgStrainEigenvector - cgStrainEigenvalue*cgStrainEigenvector:\n')
fprintf('\tmin 1: %g\t 2: %g\n',min(EigError(:,1)),min(EigError(:,2)))
fprintf('\tmax 1: %g\t 2: %g\n',max(EigError(:,1)),max(EigError(:,2)))
fprintf('\tmean 1: %g\t 2: %g\n',mean(EigError(~nanIdx,1)),mean(EigError(~nanIdx,2)))
fprintf('\tmedian 1: %g\t 2: %g\n',median(EigError(~nanIdx,1)),median(EigError(~nanIdx,2)))

%% Product of eigenvalues
prodCgStrainD = prod(cgStrainEigenvalue,2);

o(1) = 1 - max(cgStrainEigenvalue(:,1));
detCgStrain = arrayfun(@(idx)det(cgStrain(:,:,idx)),1:size(cgStrain,3));
o(2) = mean(abs(detCgStrain(~nanIdx)-1));

fprintf('lambda_1*lambda2:\n')
fprintf('\tmin = %g\n',min(prodCgStrainD))
fprintf('\tmax = %g\n',max(prodCgStrainD))
fprintf('\tmean = %g\n',mean(prodCgStrainD(~nanIdx)))
fprintf('\tmedian = %g\n',median(prodCgStrainD(~nanIdx)))
fprintf('1 - max(lambda_1) = %g\n',o(1))
fprintf('mean(abs(detCgStrain-1)) = %g\n',o(2))

end
