function [M, genes_found, gene_idx] = project_genes(genes, genes_all, pc_imputed, U)
% project_genes -- obtain gene values from compressed imputed data
%   [M, genes_found, gene_idx] = project_genes(genes, genes_all, pc_imputed, U) computes
%   gene values (M) for given gene names (genes) given all gene names (genes_all), loadings
%   (U), and imputed principal components (pc_imputed).
%
%   Since pc_imputed and U are both narrow matrices the imputed data can be
%   stored in a memory efficient way, without having to store the dense
%   matrix.

[gene_idx,locb] = ismember(lower(genes_all), lower(genes));
[~,sidx] = sort(locb(gene_idx));
idx = find(gene_idx);
idx = idx(sidx);
M = pc_imputed * U(idx,:)'; % project
genes_found = genes_all(idx);
gene_idx = find(gene_idx);
