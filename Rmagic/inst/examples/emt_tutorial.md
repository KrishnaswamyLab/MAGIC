Rmagic EMT Tutorial
================

true

<!-- README.md is generated from README.Rmd. Please edit that file -->

## MAGIC (Markov Affinity-Based Graph Imputation of Cells)

  - MAGIC imputes missing data values on sparse data sets, restoring the
    structure of the data
  - It also proves dimensionality reduction and gene expression
    visualizations
  - MAGIC can be performed on a variety of datasets
  - Here, we show the effectiveness of MAGIC on
    epithelial-to-mesenchymal transition (EMT) data

Markov Affinity-based Graph Imputation of Cells (MAGIC) is an algorithm
for denoising and transcript recover of single cells applied to
single-cell RNA sequencing data, as described in Van Dijk D *et al.*
(2018), *Recovering Gene Interactions from Single-Cell Data Using Data
Diffusion*, Cell
<https://www.cell.com/cell/abstract/S0092-8674(18)30724-4>.

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
if (!require(readr)) install.packages("readr")
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
```

    ## Loading required package: Matrix

``` r
library(readr)
library(ggplot2)
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(phateR)
```

    ## 
    ## Attaching package: 'phateR'

    ## The following object is masked from 'package:Rmagic':
    ## 
    ##     library.size.normalize

### Loading data

The example data is located in the MAGIC Github repository.

``` r
# load data
data <- read_csv("../../../data/HMLE_TGFb_day_8_10.csv.gz")
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_integer()
    ## )

    ## See spec(...) for full column specifications.

``` r
data[1:5,1:10]
```

    ## # A tibble: 5 x 10
    ##   `5S_rRNA` `5_8S_rRNA`  A1BG `A1BG-AS1`   A2M `A2M-AS1` A2ML1 `A2ML1-AS1`
    ##       <int>       <int> <int>      <int> <int>     <int> <int>       <int>
    ## 1         0           0     0          0     0         0     0           0
    ## 2         0           0     0          0     0         0     0           0
    ## 3         0           0     0          0     0         0     0           0
    ## 4         0           0     0          0     0         0     0           0
    ## 5         0           0     0          0     0         0     0           0
    ## # ... with 2 more variables: A4GALT <int>, AAAS <int>

First, we need to remove lowly expressed genes.

``` r
# keep genes expressed in at least 10 cells
keep_cols <- colSums(data > 0) > 10
data <- data[,keep_cols]
```

Ordinarily, we would remove cells with small library sizes. In this
dataset, it has already been done; however, if you wanted to do that,
you could do it with the code below.

``` r
# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(data)), bins=50) +
  geom_vline(xintercept = 1000, color='red')
```

![](emt_tutorial_files/figure-gfm/libsize_histogram-1.png)<!-- -->

``` r
if (FALSE) {
  # keep cells with at least 1000 UMIs and at most 15000
  keep_rows <- rowSums(data) > 1000 & rowSums(data) < 15000
  data <- data[keep_rows,]
}
```

We should library size normalize the data prior to MAGIC. Often we also
transform the data with either log or square root. The log transform is
commonly used, which requires adding a “pseudocount” to avoid log(0). We
normally square root instead, which has a similar form but doesn’t
suffer from instabilities at zero. For this dataset, though, it is not
necessary as the distribution of gene expression is not too extreme.

``` r
data <- library.size.normalize(data)
if (FALSE) {
  data <- sqrt(data)
}
```

### Running MAGIC

Running MAGIC is as simple as running the `magic` function. Because this
dataset is rather large, we can increase `k` from the default of 10 up
to 15.

``` r
# run MAGIC
data_MAGIC <- magic(data, k=15, genes=c("VIM", "CDH1", "ZEB1"))
```

### Results

We can plot the data before and after MAGIC to visualize the results.

``` r
ggplot(data) +
  geom_point(aes(VIM, CDH1, colour=ZEB1)) +
  scale_colour_viridis(option="B")
```

![](emt_tutorial_files/figure-gfm/plot_raw-1.png)<!-- -->

``` r
ggsave('EMT_data_R_before_magic.png', width=5, height=5)
```

``` r
ggplot(data_MAGIC) +
  geom_point(aes(VIM, CDH1, colour=ZEB1)) +
  scale_colour_viridis(option="B")
```

![](emt_tutorial_files/figure-gfm/plot_magic-1.png)<!-- -->

``` r
ggsave('EMT_data_R_after_magic.png', width=5, height=5)
```

As you can see, the gene-gene relationships are much clearer after
MAGIC.

### Visualizing MAGIC values on PCA

We can visualize the results of MAGIC on PCA with `genes="pca_only"`.

``` r
data_MAGIC_PCA <- magic(data, genes="pca_only", 
                        k=15, init=data_MAGIC)
ggplot(data_MAGIC_PCA) +
  geom_point(aes(x=PC1, y=PC2, color=data_MAGIC$result$VIM)) +
  scale_color_viridis(option="B") +
  labs(color="VIM")
```

![](emt_tutorial_files/figure-gfm/run_pca-1.png)<!-- -->

``` r
ggsave('EMT_data_R_pca_colored_by_magic.png', width=5, height=5)
```

### Using MAGIC for downstream analysis

We can look at the entire smoothed matrix with `genes='all_genes'`,
passing the original result to the argument `init` to avoid recomputing
intermediate steps. Note that this matrix may be large and could take up
a lot of memory.

``` r
data_MAGIC <- magic(data, genes="all_genes", 
                    k=15, init=data_MAGIC)
as.data.frame(data_MAGIC)[1:5, 1:10]
```

    ##          A1BG   A1BG-AS1       A2ML1      A4GALT      AAAS      AACS
    ## 1 0.001839592 0.02876726 0.008686110 0.012471223 0.1274074 0.2205920
    ## 2 0.001104917 0.03285876 0.010440814 0.006188766 0.1265774 0.2170516
    ## 3 0.002082237 0.02697853 0.008503791 0.013469548 0.1276566 0.2208905
    ## 4 0.003121141 0.03391476 0.019349338 0.017437881 0.1053139 0.2323029
    ## 5 0.001686928 0.02783747 0.007459623 0.010616114 0.1303800 0.2198455
    ##        AADAT     AAED1     AAGAB      AAK1
    ## 1 0.03046543 0.1442079 0.2745806 0.6525354
    ## 2 0.02730860 0.1485649 0.2754257 0.6392821
    ## 3 0.03098135 0.1446796 0.2742247 0.6526963
    ## 4 0.02378804 0.1464825 0.2638087 0.6597285
    ## 5 0.03099255 0.1444643 0.2769162 0.6451414

## Help

If you have any questions or require assistance using MAGIC, please
contact us at <https://krishnaswamylab.org/get-help>.
