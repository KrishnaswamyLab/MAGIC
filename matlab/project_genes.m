function [M, genes_found, lia] = project_genes(genes, pc, U, genes_all)

[lia,locb] = ismember(lower(genes_all), lower(genes));
[~,sidx] = sort(locb(lia));
idx = find(lia);
idx = idx(sidx);
%M = pc * U(idx,:)' + mu(idx);
M = pc * U(idx,:)';
genes_found = genes_all(idx);
lia = find(lia);
