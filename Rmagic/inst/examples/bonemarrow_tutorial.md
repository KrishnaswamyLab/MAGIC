Rmagic Bone Marrow Tutorial
================

true

<!-- README.md is generated from README.Rmd. Please edit that file -->

## MAGIC (Markov Affinity-Based Graph Imputation of Cells)

  - MAGIC imputes missing data values on sparse data sets, restoring the
    structure of the data
  - It also proves dimensionality reduction and gene expression
    visualizations
  - MAGIC can be performed on a variety of datasets
  - Here, we show the effectiveness of MAGIC on erythroid and myeloid
    cells developing in mouse bone marrow.

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
library(ggplot2)
library(readr)
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

In this tutorial, we will analyse myeloid and erythroid cells in mouse
bone marrow, as described in Paul et al., 2015. The example data is 
located in the PHATE Github repository and we can load it directly from 
the web.

``` r
# load data
bmmsc <- read_csv("https://github.com/KrishnaswamyLab/PHATE/raw/master/data/BMMC_myeloid.csv.gz")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_integer(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
bmmsc <- bmmsc[,2:ncol(bmmsc)]
bmmsc[1:5,1:10]
```

    ## # A tibble: 5 x 10
    ##   `0610007C21Rik;Apr… `0610007L01Rik` `0610007P08Rik;Rad2… `0610007P14Rik`
    ##                 <int>           <int>                <int>           <int>
    ## 1                   0               0                    0               0
    ## 2                   0               0                    0               1
    ## 3                   0               1                    0               2
    ## 4                   0               1                    0               1
    ## 5                   0               0                    1               0
    ## # ... with 6 more variables: `0610007P22Rik` <int>, `0610008F07Rik` <int>,
    ## #   `0610009B22Rik` <int>, `0610009D07Rik` <int>, `0610009O20Rik` <int>,
    ## #   `0610010B08Rik;Gm14434;Gm14308` <int>

First, we need to remove lowly expressed genes and cells with small
library size.

``` r
# keep genes expressed in at least 10 cells
keep_cols <- colSums(bmmsc > 0) > 10
bmmsc <- bmmsc[,keep_cols]
# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(bmmsc)), bins=50) +
  geom_vline(xintercept = 1000, color='red')
```

![](bonemarrow_tutorial_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# keep cells with at least 1000 UMIs
keep_rows <- rowSums(bmmsc) > 1000
bmmsc <- bmmsc[keep_rows,]
```

We should library size normalize and transform the data prior to MAGIC.
Many people use a log transform, which requires adding a “pseudocount”
to avoid log(0). We square root instead, which has a similar form but
doesn’t suffer from instabilities at zero.

``` r
bmmsc <- library.size.normalize(bmmsc)
bmmsc <- sqrt(bmmsc)
```

### Running MAGIC

Running MAGIC is as simple as running the `magic` function.

``` r
# run MAGIC
bmmsc_MAGIC <- magic(bmmsc, genes=c("Mpo", "Klf1", "Ifitm1"))
```

We can plot the data before and after MAGIC to visualize the results.

``` r
ggplot(bmmsc) +
  geom_point(aes(Mpo, Klf1, colour=Ifitm1)) +
  scale_colour_viridis(option="B")
```

![](bonemarrow_tutorial_files/figure-gfm/plot_raw-1.png)<!-- -->

``` r
ggsave('BMMSC_data_R_before_magic.png', width=5, height=5)
```

The data suffers from dropout to the point that we cannot infer anything
about the gene-gene relationships.

``` r
ggplot(bmmsc_MAGIC) +
  geom_point(aes(Mpo, Klf1, colour=Ifitm1)) +
  scale_colour_viridis(option="B")
```

![](bonemarrow_tutorial_files/figure-gfm/plot_magic-1.png)<!-- -->

As you can see, the gene-gene relationships are much clearer after
MAGIC. These relationships also match the biological progression we
expect to see - Ifitm1 is a stem cell marker, Klf1 is an erythroid
marker, and Mpo is a myeloid marker.

The data is a little too smooth - we can decrease `t` from the automatic
value to reduce the amount of diffusion. We pass the original result to
the argument `init` to avoid recomputing intermediate
steps.

``` r
bmmsc_MAGIC <- magic(bmmsc, genes=c("Mpo", "Klf1", "Ifitm1"), t=4, init=bmmsc_MAGIC)
ggplot(bmmsc_MAGIC) +
  geom_point(aes(Mpo, Klf1, colour=Ifitm1)) +
  scale_colour_viridis(option="B")
```

![](bonemarrow_tutorial_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave('BMMSC_data_R_after_magic.png', width=5, height=5)
```

We can look at the entire smoothed matrix with `genes='all_genes'`,
passing the original result to the argument `init` to avoid recomputing
intermediate steps. Note that this matrix may be large and could take up
a lot of memory.

``` r
bmmsc_MAGIC <- magic(bmmsc, genes="all_genes", t=4, init=bmmsc_MAGIC)
as.data.frame(bmmsc_MAGIC)[1:5, 1:10]
```

    ##   0610007C21Rik;Apr3 0610007L01Rik 0610007P08Rik;Rad26l 0610007P14Rik
    ## 1          0.1518517     0.4808969           0.05179379     0.3415315
    ## 2          0.1230258     0.3904692           0.12708682     0.3880762
    ## 3          0.1344016     0.5496261           0.06617484     0.3745274
    ## 4          0.1132422     0.3812295           0.11713295     0.3683106
    ## 5          0.1417958     0.4637729           0.04572850     0.3465252
    ##   0610007P22Rik 0610009B22Rik 0610009D07Rik 0610009O20Rik
    ## 1    0.02452145    0.05066031    0.06728059     0.1611080
    ## 2    0.04026323    0.07548813    0.11778619     0.3763095
    ## 3    0.02323954    0.05561339    0.05379475     0.1802045
    ## 4    0.03597300    0.06506830    0.08923059     0.3736334
    ## 5    0.02303185    0.05190157    0.05275807     0.1621462
    ##   0610010F05Rik;mKIAA1841;Kiaa1841 0610010K14Rik;Rnasek
    ## 1                       0.03023858            1.0447754
    ## 2                       0.03823837            0.9875507
    ## 3                       0.04074101            1.0408079
    ## 4                       0.03201622            1.0149211
    ## 5                       0.02737464            0.9940195

### Visualizing MAGIC values on PCA

We can visualize the results of MAGIC on PCA as follows.

``` r
bmmsc_MAGIC_PCA <- as.data.frame(prcomp(bmmsc_MAGIC)$x)
ggplot(bmmsc_MAGIC_PCA) +
  geom_point(aes(x=PC1, y=PC2, color=bmmsc_MAGIC$result$Klf1)) +
  scale_color_viridis(option="B") +
  labs(color="Klf1")
```

![](bonemarrow_tutorial_files/figure-gfm/run_pca-1.png)<!-- -->

``` r
ggsave('BMMSC_data_R_pca_colored_by_magic.png', width=5, height=5)
```

### Visualizing MAGIC values on PHATE

We can visualize the results of MAGIC on PHATE as follows.

``` r
bmmsc_PHATE <- phate(bmmsc, k=4, a=100, t=20)
ggplot(bmmsc_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=bmmsc_MAGIC$result$Klf1)) +
  scale_color_viridis(option="B") +
  labs(color="Klf1")
```

![](bonemarrow_tutorial_files/figure-gfm/run_phate-1.png)<!-- -->

``` r
ggsave('BMMSC_data_R_phate_colored_by_magic.png', width=5, height=5)
```
