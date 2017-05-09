function data_imputed = run_magic(data, t, varargin)
% MAGIC (Markov Affinity-based Graph Imputation of Cells)
% data must have cells on the rows and genes on the columns
% t is diffusion time
% varargin:
%   'npca' (default = 20)
%       perform fast random PCA before computing distances
%   'ka' (default = 10)
%       k of adaptive kernel
%       0 for non-adaptive (standard gaussian) kernel with bandwidth
%       epsilon
%   'k' (default = 30)
%       k of kNN graph
%   'epsilon' (default = 1)
%       kernel bandwith, if epsilon = 0 kernel will be uniform, i.e.
%       unweighted kNN graph (ka will be ignored)
%   'rescale_to' (default = 0, no rescale)
%       rescale genes to 'rescale_to' percentile
%       set to 0 for log scaled data

% set up default parameters
k = 30;
ka = 10;
npca = 20;
rescale_to = 0;
epsilon = 1;

% get the input parameters
if ~isempty(varargin)
    for j = 1:length(varargin)
        % k nearest neighbor 
        if strcmp(varargin{j}, 'ka')
            ka = varargin{j+1};
        end
        % for knn-autotune
        if strcmp(varargin{j}, 'k')
            k = varargin{j+1};
        end
        % epsilon
        if strcmp(varargin{j}, 'epsilon')
            epsilon = varargin{j+1};
        end
        % npca to project data
        if strcmp(varargin{j}, 'npca')
            npca = varargin{j+1};
        end
        % sigma of kernel bandwidth
        if strcmp(varargin{j}, 'rescale_to')
            rescale_to = varargin{j+1};
        end
    end
end

N = size(data, 1); % number of cells

disp 'PCA'
data_centr = bsxfun(@minus, data, mean(data,1));
[U,~,~] = randPCA(data_centr', npca); % fast random svd
%[U,~,~] = svds(data', npca);
data_pc = data_centr * U; % PCA project

disp 'Computing distances'
[idx, dist] = knnsearch(data_pc, data_pc, 'k', k);

disp 'Adapting sigma'
dist = bsxfun(@rdivide, dist, dist(:,ka));

i = repmat((1:N)',1,size(idx,2));
i = i(:);
j = idx(:);
s = dist(:);
if epsilon > 0
    W = sparse(i, j, s);
else
    W = sparse(i, j, ones(size(s))); % unweighted kNN graph
end

disp 'Symmetrize distances'
W = W + W';

if epsilon > 0
    disp 'Computing kernel'
    [i,j,s] = find(W);
    i = [i; (1:N)'];
    j = [j; (1:N)'];
    s = [s./(epsilon^2); zeros(N,1)];
    s = exp(-s);
    W = sparse(i,j,s);
end

disp 'Markov normalization'
W = bsxfun(@rdivide, W, sum(W,2)); % Markov normalization

W = full(W);

disp(['Diffusing for ' num2str(t) ' steps']);
W_t = W^t; % diffuse

disp 'Imputing'
data_imputed = W_t * data; % impute

% Rescale
if rescale_to > 0
    if ~any(data(:)<0)
        disp 'Rescaling'
        MR = prctile(data, rescale_to);
        M = max(data);
        MR(MR == 0) = M(MR == 0);
        MR_new = prctile(data_imputed, rescale_to);
        M_new = max(data_imputed);
        MR_new(MR_new == 0) = M_new(MR_new == 0);
        max_ratio = MR ./ MR_new;
        data_imputed = data_imputed .* repmat(max_ratio, N, 1);
    else
        disp('Negative values detected (log scaled data?) so no rescale is done.')
    end
end

disp 'done'
