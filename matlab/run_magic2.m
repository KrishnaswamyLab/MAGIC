function data_imputed = run_magic2(data, t, varargin)
% run MAGIC
%  Implemetation that does not use Mauro Maggioni's Diffusion Geometry package

% data must have cells on the rows and genes on the columns
% t is diffusion time
% varargin:
%   'npca' (defauklt = 20)
%       perform fast random PCA before computing distances
%   'ka' (default = 10)
%       k of adaptive kernel
%       0 for non-adaptive (standard gaussian) kernel with bandwidth sigma
%   'k' (defauklt = 30)
%       k of kNN graph
%   'rescale_to' (default = 0, no rescale)
%       rescale genes to 'rescale_to' percentile
%       set to 0 for log scaled data

% set up default parameters
k = 30;
ka = 10;
npca = 20;
rescale_to = 0;

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

disp 'PCA'
[U,~,~] = randPCA(data', npca); % fast random svd
%[U,~,~] = svds(data', npca);
data_pc = data * U; % PCA project

disp 'Computing kernel'
D = squareform(pdist(data_pc)); % compute distances
knnD = sort(D); % get kNN
th = knnD(k+1,:); % find knn cutoff
ind_zero = D > repmat(th',1,size(D,2)); % kNN cutoff

disp 'Adapting sigma'
sigma = knnD(ka+1,:); % locally adaptive sigma
D = bsxfun(@rdivide,D,sigma); % adapt distances
W = exp(-D.^2); % gaussian kernel
W = W + W'; % symmetrize
W(ind_zero) = 0; % kNN cutoff
L = bsxfun(@rdivide, W, sum(W,2)); % Markov normalization

disp(['Diffusing for ' num2str(t) ' steps']);
L_t = L^t; % diffuse

disp 'Imputing'
data_imputed = L_t * data; % impute

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
        data_imputed = data_imputed .* repmat(max_ratio, size(data,1), 1);
    else
        disp('Negative values detected (log scaled?) so no rescale is done.')
    end
end

disp 'done'
