function [t_opt, r2_vec] = compute_optimal_t(data, DiffOp, varargin)
% [t_opt, r2_vec] = compute_optimal_t(data, DiffOp, varargin)
%   data - input data
%   DiffOp - diffusion operator
%   varargin:
%       t_max - max t to try
%       n_genes - number of random genes to compute optimal t on, should be
%        at least 100, fewer is faster
%       make_plots - draw convergence as a function of t with which we
%        select the optimal t

t_max = 32;
n_genes = size(data,2);
make_plots = true;

if ~isempty(varargin)
    for j = 1:length(varargin)
        if strcmp(varargin{j}, 't_max')
            t_max = varargin{j+1};
        end
        if strcmp(varargin{j}, 'n_genes')
            n_genes = varargin{j+1};
        end
        if strcmp(varargin{j}, 'make_plots')
            make_plots = varargin{j+1};
        end
    end
end

if ~issparse(DiffOp)
    DiffOp = sparse(DiffOp);
end

if n_genes > size(data,2)
    disp 'n_genes too large, capping n_genes at maximum possible number of genes'
    n_genes = size(data,2)
end

idx_genes = randsample(size(data,2), n_genes);
data_imputed = data;
data_imputed = data_imputed(:,idx_genes);

if min(data_imputed(:)) < 0
    disp 'data has negative values, shifting to positive'
    data_imputed = data_imputed - min(data_imputed(:));
end

r2_vec = nan(t_max,1);
data_prev = data_imputed;
data_prev = bsxfun(@rdivide, data_prev, sum(data_prev));
disp 'computing optimal t'
for I=1:t_max
    data_imputed = DiffOp * data_imputed;
    data_curr = data_imputed;
    data_curr = bsxfun(@rdivide, data_curr, sum(data_curr));
    r2 = rsquare(data_prev(:), data_curr(:));
    r2_vec(I) = 1 - r2;
    data_prev = data_curr;
end

t_opt = find(r2_vec < 0.05, 1, 'first') + 1;

disp(['optimal t = ' num2str(t_opt)]);

if make_plots
    figure;
    hold all;
    plot(1:t_max, r2_vec, '*-');
    plot(t_opt, r2_vec(t_opt), 'or', 'markersize', 10);
    xlabel 't'
    ylabel '1 - R^2(data_{t},data_{t-1})'
    axis tight
    ylim([0 1]);
    plot(xlim, [0.05 0.05], '--k');
    legend({'y' 'optimal t' 'y=0.05'});
    set(gca,'xtick',1:t_max);
    set(gca,'ytick',0:0.1:1);
end


