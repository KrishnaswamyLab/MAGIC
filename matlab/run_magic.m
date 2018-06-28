function [pc_t, U] = run_magic(data, varargin)

npca = 100;
k = 15;
a = 15;
t = [];
distfun = 'euclidean';

% get input parameters
for i=1:length(varargin)
    % k for knn adaptive sigma
    if(strcmp(varargin{i},'k'))
       k = lower(varargin{i+1});
    end
    % a (alpha) for alpha decaying kernel
    if(strcmp(varargin{i},'a'))
       a = lower(varargin{i+1});
    end
    % diffusion time
    if(strcmp(varargin{i},'t'))
       t = lower(varargin{i+1});
    end
    % npca
    if(strcmp(varargin{i},'npca'))
       npca = lower(varargin{i+1});
    end
end

% PCA
disp 'doing PCA'
[U,~,~] = randPCA(data', npca);
pc = data * U;

% compute kernel
disp 'computing kernel'
K = compute_kernel(pc, 'k', k, 'a', a, 'distfun', distfun);

% row stochastic
P = bsxfun(@rdivide, K, sum(K,2));

% optimal t
if isempty(t)
    t = compute_optimal_t(pc, P);
end

% MAGIC on PCA
disp 'imputing'
pc_t = pc;
for I=1:t
    disp(['t = ' num2str(I)]);
    pc_t = P * pc_t;
end

disp 'done.'
