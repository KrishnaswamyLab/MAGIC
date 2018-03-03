%% load data (e.g. 10x data)
% data should be cells as rows and genes as columns
%sample_dir = 'path_to_data/';
%[data, gene_names, gene_ids, cells] = load_10x(sample_dir);

%% load EMT data
file = 'EMT.csv'; % unzip EMT.csv.zip
data = importdata(file);
gene_names = data.colheaders;
data = data.data;

%% MAGIC
npca = 20; % ususally between 10 and 200
ka = 10; % can be smaller, eg 3
k = 30; % can be smaller, eg 9 (3*ka)
t = []; % automatically find and use optimal t
lib_size_norm = true; % library size normalize
log_transform = false; % log transform, some data requires this
[data_imputed, DiffOp] = run_magic(data, t, 'npca', npca, 'ka', ka, 'k', k, ...
    'lib_size_norm', lib_size_norm, 'log_transform', log_transform);

%% plot
plot_genes = {'Cdh1', 'Vim', 'Fn1', 'Zeb1'};
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
x = data_imputed(:, ismember(lower(gene_names), lower(plot_genes{1})));
y = data_imputed(:, ismember(lower(gene_names), lower(plot_genes{2})));
z = data_imputed(:, ismember(lower(gene_names), lower(plot_genes{3})));
c = data_imputed(:, ismember(lower(gene_names), lower(plot_genes{4})));
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




