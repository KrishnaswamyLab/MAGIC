function [pc_imputed, U, pc] = run_magic(data, varargin)
% run_magic  Run MAGIC for imputing and denoising of single-cell data
%   [pc_imputed, U, pc] = run_magic(data, varargin) runs MAGIC on data (rows:
%   cells, columns: genes) with default parameter settings and returns the
%   imputed data in a compressed format.
%
%   The compressed format consists of loadings (U) and imputed principal
%   components (pc_imputed). To obtain gene values form this compressedf format
%   either run project_genes.m or manually project (pc_imputed * U') either all
%   genes or a subset (pc_imputed * U(idx,:)'). Also returned are the original
%   principal components (pc);
%
%   Since pc_imputed and U are both narrow matrices the imputed data can be
%   stored in a memory efficient way, without having to store the dense
%   matrix.
%
%   Supplied data can be a sparse matrix, in which case MAGIC will be more
%   memory efficient.
%
%   [...] = phate(data, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to
%   specify optional parameter name/value pairs that control further details
%   of PHATE.  Parameters are:
%
%   'npca' - number of PCA components to do MAGIC on. Defaults to 100.
%
%   'k' - number of nearest neighbors for bandwidth of adaptive alpha
%   decaying kernel or, when a=[], number of nearest neighbors of the knn
%   graph. For the unweighted kernel we recommend k to be a bit larger,
%   e.g. 10 or 15. Defaults to 10.
%
%   'a' - alpha of alpha decaying kernel. when a=[] knn (unweighted) kernel
%   is used. Defaults to 15.
%
%   't' - number of diffusion steps. Defaults to [] wich autmatically picks
%   the optimal t.
%
%   'distfun' - Distance function used to compute kernel. Defaults to
%   'euclidean'.
%
%   'make_plot_opt_t' - Boolean flag for plotting the optimal t analysis.
%   Defaults to true.

npca = 100;
k = 10;
a = 15;
t = [];
distfun = 'euclidean';
make_plot_opt_t = true;

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
    % make plot optimal t
    if(strcmp(varargin{i},'make_plot_opt_t'))
       make_plot_opt_t = lower(varargin{i+1});
    end
end

% PCA
disp 'doing PCA'
[U,~,~] = randPCA(data', npca); % this is svd
pc = data * U; % this is PCA without mean centering to be able to handle sparse data

% compute kernel
disp 'computing kernel'
K = compute_kernel(pc, 'k', k, 'a', a, 'distfun', distfun);

% row stochastic
P = bsxfun(@rdivide, K, sum(K,2));

% optimal t
if isempty(t)
    disp 'imputing using optimal t'
    pc_imputed = compute_optimal_t(pc, P, 'make_plot', make_plot_opt_t);
else
    disp 'imputing using provided t'
    pc_imputed = pc;
    for I=1:t
        disp(['t = ' num2str(I)]);
        pc_imputed = P * pc_imputed;
    end
end

disp 'done.'
