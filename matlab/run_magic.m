function data_imputed = run_magic(data, t, varargin)
% run MAGIC

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
%   'sigma' (default = 1)
%       sigma of kernel bandwidth
%       for adaptive kernel (ka>0) the effective bandwidth is sigma * the
%       distance to the ka-th neighbor
%   'rescale_to' (default = 0, no rescale)
%       rescale genes to 'rescale_to' percentile
%       set to 0 for log scaled data

% set up default parameters
k = 30;
ka = 10;
npca = 20;
sigma = 1;
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
        if strcmp(varargin{j}, 'sigma')
            sigma = varargin{j+1};
        end
        % sigma of kernel bandwidth
        if strcmp(varargin{j}, 'rescale_to')
            rescale_to = varargin{j+1};
        end
    end
end

% Kernel
disp 'Computing kernel'
Options.Display = 1;
Options.Epsilon = sigma;
Options.kNN = k;
Options.kNNAutotune = ka;
Options.NNMaxDim = npca;
Options.Normalization = 'markov';
G = GraphDiffusion(data', 0, Options);
L = full(G.T);

% Diffuse
disp(['Diffusing for ' num2str(t) ' steps']);
L_t = L^t;

% Impute
disp 'Imputing'
data_imputed = L_t * data;

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
