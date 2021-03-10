function Y = svdpca(X, k, method)

if ~exist('method','var')
    method = 'svd';
end

X = bsxfun(@minus, X, mean(X));

switch lower(method)
    case 'svd'
        disp 'PCA using SVD'
        [U,~,~] = svds(X', k);
        Y = X * U;
    case 'random'
        disp 'PCA using random SVD'
        [U,~,~] = randPCA(X', k);
        Y = X * U;
    case 'lra'
        disp 'LRA using random SVD'
        [U,S,V] = svds(X, k);
        Y = U*S*V';
    case 'lra_random'
        disp 'LRA using random SVD'
        [U,S,V] = randPCA(X, k);
        Y = U*S*V';
end
