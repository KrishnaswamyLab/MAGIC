function [data, gene_names, gene_ids, cells] = load_10x(data_dir, varargin)
% [data, gene_names, gene_ids, cells] = load_10x(data_dir, varargin)
%   loads 10x sparse format data
%   data_dir is dir that contains matrix.mtx, genes.tsv and barcodes.tsv
%   varargin
%       'sparse', true -- returns data matrix in sparse format (default 'false')

return_sparse = false;

if isempty(data_dir)
    data_dir = './';
elseif data_dir(end) ~= '/'
    data_dir = [data_dir '/'];
end

for i=1:length(varargin)-1
    if (strcmp(varargin{i}, 'sparse'))
        return_sparse = varargin{i+1};
    end
end

filename_dataMatrix = [data_dir 'matrix.mtx'];
filename_genes = [data_dir 'genes.tsv'];
filename_cells = [data_dir 'barcodes.tsv'];


% Read in gene expression matrix (sparse matrix)
% Rows = genes, columns = cells
fprintf('LOADING\n')
dataMatrix = mmread(filename_dataMatrix);
fprintf('  Data matrix (%i cells x %i genes): %s\n', ...
        size(dataMatrix'), ['''' filename_dataMatrix '''' ])

% Read in row names (gene names / IDs)
dataMatrix_genes = table2cell( ...
                   readtable(filename_genes, ...
                             'FileType','text','ReadVariableNames',0));
dataMatrix_cells = table2cell( ...
                   readtable(filename_cells, ...
                             'FileType','text','ReadVariableNames',0));

% Remove empty cells
col_keep = any(dataMatrix,1);
dataMatrix = dataMatrix(:,col_keep);
dataMatrix_cells = dataMatrix_cells(col_keep,:);
fprintf('  Removed %i empty cells\n', full(sum(~col_keep)))

% Remove empty genes
genes_keep = any(dataMatrix,2);
dataMatrix = dataMatrix(genes_keep,:);
dataMatrix_genes = dataMatrix_genes(genes_keep,:);
fprintf('  Removed %i empty genes\n', full(sum(~genes_keep)))

data = dataMatrix';
if ~return_sparse
    data = full(data);
end
gene_names = dataMatrix_genes(:,2);
gene_ids = dataMatrix_genes(:,1);
cells = dataMatrix_cells;
