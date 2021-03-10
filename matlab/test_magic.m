%% load data (e.g. 10x data)
% data should be cells as rows and genes as columns
%sample_dir = 'path_to_data/';
%[data, gene_names, gene_ids, cells] = load_10x(sample_dir);

%% load EMT data
file = '../data/HMLE_TGFb_day_8_10.csv'; %% gunzip ../data/HMLE_TGFb_day_8_10.csv.gz
data = importdata(file);
gene_names = data.colheaders;
data = data.data;

%% library size normalization
libsize = sum(data,2);
data = bsxfun(@rdivide, data, libsize) * median(libsize);

%% log transform -- usually one would log transform the data. Here we don't do it.
%data = log(data + 0.1);

%% MAGIC
[pc_imputed, U, pc] = run_magic(data, 'npca', 100, 'k', 15, 'a', 15, 'make_plot_opt_t', true);

%% project genes
plot_genes = {'Cdh1', 'Vim', 'Fn1', 'Zeb1'};
[M_imputed, genes_found] = project_genes(plot_genes, gene_names, pc_imputed, U);

%% plot
ms = 20;
v = [-45 20];
% before MAGIC
x = data(:, ismember(lower(gene_names), lower(plot_genes{1})));
y = data(:, ismember(lower(gene_names), lower(plot_genes{2})));
z = data(:, ismember(lower(gene_names), lower(plot_genes{3})));
c = data(:, ismember(lower(gene_names), lower(plot_genes{4})));
figure;
subplot(2,2,1);
scatter(y, x, ms, c, 'filled');
colormap(parula);
axis tight
xlabel(plot_genes{2});
ylabel(plot_genes{1});
h = colorbar;
%ylabel(h,plot_genes{4});
title 'Before MAGIC'

subplot(2,2,2);
scatter3(x, y, z, ms, c, 'filled');
colormap(parula);
axis tight
xlabel(plot_genes{1});
ylabel(plot_genes{2});
zlabel(plot_genes{3});
%h = colorbar;
ylabel(h,plot_genes{4});
view(v);
title 'Before MAGIC'

% plot after MAGIC
x = M_imputed(:,1);
y = M_imputed(:,2);
z = M_imputed(:,3);
c = M_imputed(:,4);
subplot(2,2,3);
scatter(y, x, ms, c, 'filled');
colormap(parula);
axis tight
xlabel(plot_genes{2});
ylabel(plot_genes{1});
h = colorbar;
%ylabel(h,plot_genes{4});
title 'After MAGIC'

subplot(2,2,4);
scatter3(x, y, z, ms, c, 'filled');
colormap(parula);
axis tight
xlabel(plot_genes{1});
ylabel(plot_genes{2});
zlabel(plot_genes{3});
%h = colorbar;
ylabel(h,plot_genes{4});
view(v);
title 'After MAGIC'

%% plot PCA before MAGIC
figure;
c = data(:, ismember(lower(gene_names), lower(plot_genes{4})));
Y = svdpca(pc, 3, 'random'); % original PCs are not mean centered so doing proper PCA here
% alternative is to do proper PCA on data:
%Y = svdpca(data, 3, 'random');
scatter3(Y(:,1), Y(:,2), Y(:,3), ms, c, 'filled');
colormap(parula);
axis tight
xlabel 'PC1'
ylabel 'PC2'
zlabel 'PC3'
h = colorbar;
ylabel(h,plot_genes{4});
view([-50 22]);
title 'Before MAGIC'

%% plot PCA after MAGIC
figure;
c = M_imputed(:,4);
Y = svdpca(pc_imputed, 3, 'random'); % original PCs are not mean centered so doing proper PCA here
% alternative is to go to full imputed data and then do proper PCA:
%data_imputed = pc_imputed * U'; % project full data
%Y = svdpca(data_imputed, 3, 'random');
scatter3(Y(:,1), Y(:,2), Y(:,3), ms, c, 'filled');
colormap(parula);
axis tight
xlabel 'PC1'
ylabel 'PC2'
zlabel 'PC3'
h = colorbar;
ylabel(h,plot_genes{4});
view([-50 22]);
title 'After MAGIC'
