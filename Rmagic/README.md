Rmagic v1.0.0
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis CI
Build](https://api.travis-ci.com/KrishnaswamyLab/MAGIC.svg?branch=master)](https://travis-ci.com/KrishnaswamyLab/MAGIC)
[![Read the
Docs](https://img.shields.io/readthedocs/magic.svg)](https://magic.readthedocs.io/)
[![Cell Publication
DOI](https://zenodo.org/badge/DOI/10.1016/j.cell.2018.05.061.svg)](https://doi.org/10.1016/j.cell.2018.05.061)
[![Twitter](https://img.shields.io/twitter/follow/KrishnaswamyLab.svg?style=social&label=Follow)](https://twitter.com/KrishnaswamyLab)
[![Github
Stars](https://img.shields.io/github/stars/KrishnaswamyLab/MAGIC.svg?style=social&label=Stars)](https://github.com/KrishnaswamyLab/MAGIC/)

Markov Affinity-based Graph Imputation of Cells (MAGIC) is an algorithm
for denoising and transcript recover of single cells applied to
single-cell RNA sequencing data, as described in Van Dijk D *et al.*
(2018), *Recovering Gene Interactions from Single-Cell Data Using Data
Diffusion*, Cell
<https://www.cell.com/cell/abstract/S0092-8674(18)30724-4>.

  - MAGIC imputes missing data values on sparse data sets, restoring the
    structure of the data
  - It also proves dimensionality reduction and gene expression
    visualizations
  - MAGIC can be performed on a variety of datasets
  - Here, we show the usage of MAGIC on a toy dataset
  - You can view further examples of MAGIC on real data in our notebooks
    under
        `inst/examples`:
      - <https://github.com/KrishnaswamyLab/MAGIC/blob/master/Rmagic/inst/examples/EMT_tutorial.md>
      - <https://github.com/KrishnaswamyLab/MAGIC/blob/master/Rmagic/inst/examples/bonemarrow_tutorial.md>

### Installation

To use MAGIC, you will need to install both the R and Python packages.

In R, run these commands to install MAGIC and all dependencies:

``` r
if (!require(devtools)) install.packages(devtools)
if (!require(Rmagic)) devtools::install_github("KrishnaswamyLab/magic/Rmagic")
```

In a terminal, run the following command to install the Python
repository.

``` bash
pip install --user git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python
```

We’ll install a couple more tools for this tutorial.

``` r
if (!require(viridis)) install.packages("viridis")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(phateR)) install.packages("phateR")
```

If you have never used PHATE, you should also install PHATE from the
command line as follows:

``` bash
pip install --user phate
```

### Loading packages

We load the Rmagic package and a few others for convenience functions.

``` r
library(Rmagic)
#> Loading required package: Matrix
library(ggplot2)
library(viridis)
#> Loading required package: viridisLite
library(phateR)
#>
#> Attaching package: 'phateR'
#> The following object is masked from 'package:Rmagic':
#>
#>     library.size.normalize
```

### Loading data

The example data is located in the MAGIC Github repository.

``` r
# load data
data(magic_testdata)
magic_testdata[1:5,1:10]
#>       A1BG-AS1     AAMDC      AAMP AARSD1 ABCA12 ABCG2    ABHD13
#> 6564 0.0000000 0.0000000 0.0000000      0      0     0 0.0000000
#> 3835 0.0000000 0.8714711 0.0000000      0      0     0 0.8714711
#> 6318 0.7739207 0.0000000 0.7739207      0      0     0 0.0000000
#> 3284 0.0000000 0.0000000 0.0000000      0      0     0 0.0000000
#> 1171 0.0000000 0.0000000 0.0000000      0      0     0 0.0000000
#>      AC007773.2 AC011998.4 AC013470.6
#> 6564          0          0          0
#> 3835          0          0          0
#> 6318          0          0          0
#> 3284          0          0          0
#> 1171          0          0          0
```

### Running MAGIC

Running MAGIC is as simple as running the `magic` function.

``` r
# run MAGIC
data_MAGIC <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"))
```

We can plot the data before and after MAGIC to visualize the results.

``` r
ggplot(magic_testdata) +
  geom_point(aes(VIM, CDH1, colour=ZEB1)) +
  scale_colour_viridis(option="B")
```

<img src="man/figures/README-plot_raw-1.png" width="100%" />

The data suffers from dropout to the point that we cannot infer anything
about the gene-gene relationships.

``` r
ggplot(data_MAGIC) +
  geom_point(aes(VIM, CDH1, colour=ZEB1)) +
  scale_colour_viridis(option="B")
```

<img src="man/figures/README-plot_magic-1.png" width="100%" />

As you can see, the gene-gene relationships are much clearer after
MAGIC.

The data is sometimes a little too smooth - we can decrease `t` from the
automatic value to reduce the amount of diffusion. We pass the original
result to the argument `init` to avoid recomputing intermediate
steps.

``` r
data_MAGIC <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"), t=4, init=data_MAGIC)
ggplot(data_MAGIC) +
  geom_point(aes(VIM, CDH1, colour=ZEB1)) +
  scale_colour_viridis(option="B")
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

We can look at the entire smoothed matrix with `genes='all_genes'`,
passing the original result to the argument `init` to avoid recomputing
intermediate steps. Note that this matrix may be large and could take up
a lot of
memory.

``` r
data_MAGIC <- magic(magic_testdata, genes="all_genes", t=4, init=data_MAGIC)
as.data.frame(data_MAGIC)[1:5, 1:10]
#>        A1BG-AS1      AAMDC      AAMP     AARSD1     ABCA12      ABCG2
#> 6564 0.02318932 0.06636833 0.1814952 0.01998665 0.03509270 0.01091799
#> 3835 0.02305522 0.06726339 0.1597559 0.01630487 0.02857813 0.01335037
#> 6318 0.02640804 0.06546350 0.1934657 0.01469336 0.03087314 0.01564654
#> 3284 0.02287164 0.06194795 0.1637859 0.01656443 0.02939938 0.01340856
#> 1171 0.02824691 0.06262350 0.1767155 0.01483140 0.03406291 0.01458588
#>          ABHD13   AC007773.2   AC011998.4  AC013470.6
#> 6564 0.06926381 0.0007282491 0.0012694154 0.003669113
#> 3835 0.06487983 0.0007331720 0.0016187119 0.002646376
#> 6318 0.07547549 0.0012984700 0.0020315661 0.005617329
#> 3284 0.06620499 0.0007615835 0.0007540187 0.001993851
#> 1171 0.07200313 0.0012521371 0.0025994098 0.005028384
```

### Visualizing MAGIC values on PCA

We can visualize the results of MAGIC on PCA as follows.

``` r
data_MAGIC_PCA <- as.data.frame(prcomp(data_MAGIC)$x)
ggplot(data_MAGIC_PCA) +
  geom_point(aes(x=PC1, y=PC2, color=data_MAGIC$result$VIM)) +
  scale_color_viridis(option="B") +
  labs(color="VIM")
```

<img src="man/figures/README-run_pca-1.png" width="100%" />

### Visualizing MAGIC values on PHATE

We can visualize the results of MAGIC on PHATE as follows. We set `t`
and `k` manually, because this toy dataset is really too small to make
sense with PHATE; however, the default values work well for single-cell
genomic data.

``` r
data_PHATE <- phate(magic_testdata, k=3, t=15)
ggplot(data_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=data_MAGIC$result$VIM)) +
  scale_color_viridis(option="B") +
  labs(color="VIM")
```

<img src="man/figures/README-run_phate-1.png" width="100%" />
