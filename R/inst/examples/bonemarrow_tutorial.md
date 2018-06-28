Rmagic Bone Marrow Tutorial
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## MAGIC (Markov Affinity-Based Graph Imputation of Cells)

  - MAGIC imputes missing data values on sparse data sets, restoring the
    structure of the data
  - It also proves dimensionality reduction and gene expression
    visualizations
  - MAGIC can be performed on a variety of datasets
  - Here, we show the effectiveness of MAGIC on
    epithelial-to-mesenchymal transition (EMT) data

To use MAGIC, cite: MAGIC: A diffusion-based imputation method reveals
gene-gene interactions in single-cell RNA-sequencing data Biorxiv
(preprint) February 2017. DOI: doi.org/10.1101/111591

### Installation

To use MAGIC, you will need to install both the R and Python packages.

In R, run these commands to install MAGIC and all dependencies:

``` r
if (!require(devtools)) install.packages(devtools)
if (!require(Rmagic)) devtools::install_github("KrishnaswamyLab/magic/R")
```

In a terminal, run the following command to install the Pythonn
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

    ## Warning: package 'Matrix' was built under R version 3.4.4

``` r
library(ggplot2)
library(readr)
library(viridis)
```

    ## Warning: package 'viridis' was built under R version 3.4.4

    ## Loading required package: viridisLite

    ## Warning: package 'viridisLite' was built under R version 3.4.4

``` r
library(phateR)
```

    ## 
    ## Attaching package: 'phateR'

    ## The following object is masked from 'package:Rmagic':
    ## 
    ##     library.size.normalize

### Loading data

The example data is located in the PHATE Github repository. We can load
it directly from the web.

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
    ##   `0610007C21Rik;Apr~ `0610007L01Rik` `0610007P08Rik;Rad2~ `0610007P14Rik`
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
hist(rowSums(bmmsc))
```

![](bonemarrow_tutorial_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# keep cells with at least 2000 UMIs
keep_rows <- rowSums(bmmsc) > 2000
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
  scale_colour_viridis()
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
  scale_colour_viridis()
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
bmmsc_MAGIC <- magic(bmmsc, genes=c("Mpo", "Klf1", "Ifitm1"), t=8, init=bmmsc_MAGIC)
ggplot(bmmsc_MAGIC) +
  geom_point(aes(Mpo, Klf1, colour=Ifitm1)) + 
  scale_colour_viridis()
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
bmmsc_MAGIC <- magic(bmmsc, genes="all_genes", init=bmmsc_MAGIC)
as.data.frame(bmmsc_MAGIC)[1:5, 1:10]
```

    ##   0610007C21Rik;Apr3 0610007L01Rik 0610007P08Rik;Rad26l 0610007P14Rik
    ## 1          0.1751866     0.5540753           0.06156801     0.4074682
    ## 2          0.1546742     0.4319843           0.11617584     0.4282041
    ## 3          0.1687694     0.5702529           0.06386686     0.4172185
    ## 4          0.1494420     0.4273008           0.11111974     0.4211829
    ## 5          0.1708924     0.5635437           0.06133188     0.4162455
    ##   0610007P22Rik 0610009B22Rik 0610009D07Rik 0610009O20Rik
    ## 1    0.03024783    0.05645037    0.07368982     0.1979285
    ## 2    0.03774436    0.08081385    0.12355306     0.3909738
    ## 3    0.02894124    0.05833093    0.07277939     0.2002494
    ## 4    0.03579443    0.07634111    0.11043661     0.3862620
    ## 5    0.02948439    0.05587387    0.07075943     0.2007730
    ##   0610010F05Rik;mKIAA1841;Kiaa1841 0610010K14Rik;Rnasek
    ## 1                       0.03743578             1.160120
    ## 2                       0.03603500             1.090403
    ## 3                       0.03885539             1.150197
    ## 4                       0.03315765             1.096635
    ## 5                       0.03766863             1.155895

### Visualizing MAGIC values on PHATE

We can visualize the results of MAGIC on PHATE as follows.
mKIAA1317;Kctd16

``` r
bmmsc_PHATE <- phate(bmmsc)
ggplot(bmmsc_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=bmmsc_MAGIC$result$Klf1)) +
  scale_color_viridis() +
  labs(color="Mpo")
```

![](bonemarrow_tutorial_files/figure-gfm/run_phate-1.png)<!-- -->

``` r
ggsave('BMMSC_data_R_phate_colored_by_magic.png', width=5, height=5)
```
