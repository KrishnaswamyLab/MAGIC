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

We should library size normalize and transform the data prior to MAGIC.
Many people use a log transform, which requires adding a “pseudocount”
to avoid log(0). We square root instead, which has a similar form but
doesn’t suffer from instabilities at zero.

``` r
bmmsc <- library.size.normalize(bmmsc)
bmmsc <- sqrt(bmmsc)
```

### Running MAGIC

Running MAGIC is as simple as running the `magic` function. We use
`rescale_percent=99` to make the output range comparable to the input,
but this is not strictly necessary for downstream analysis.

``` r
# run MAGIC
bmmsc_MAGIC <- magic(bmmsc, genes=c("Mpo", "Klf1", "Ifitm1"), rescale_percent=99)
```

### Results

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

``` r
ggsave('BMMSC_data_R_after_magic.png', width=5, height=5)
```

As you can see, the gene-gene relationships are much clearer after
MAGIC. These relationships also match the biological progression we
expect to see - Ifitm1 is a stem cell marker, Klf1 is an erythroid
marker, and Mpo is a myeloid marker.

We can look at the entire smoothed matrix with `genes='all_genes'`,
passing the original result to the argument `init` to avoid recomputing
intermediate steps. Note that this matrix may be large and could take up
a lot of
memory.

``` r
bmmsc_MAGIC <- magic(bmmsc, genes="all_genes", rescale_percent=99, init=bmmsc_MAGIC)
as.data.frame(bmmsc_MAGIC)[1:5, 1:10]
```

    ##   0610007C21Rik;Apr3 0610007L01Rik 0610007P08Rik;Rad26l 0610007P14Rik
    ## 1           1.860502      1.690168            1.1373480      1.905561
    ## 2           2.109082      2.176668            0.9942011      2.040323
    ## 3           2.156509      2.015408            2.0982907      2.276530
    ## 4           2.137702      2.262988            1.0380500      2.076257
    ## 5           2.066830      1.971425            1.9921411      2.203174
    ##   0610007P22Rik 0610008F07Rik 0610009B22Rik 0610009D07Rik 0610009O20Rik
    ## 1      1.262215             0      1.416935     0.9999595      1.498086
    ## 2      1.278548             0      1.295572     1.0363637      1.207875
    ## 3      1.966309             0      1.844326     1.7560429      2.438793
    ## 4      1.302293             0      1.320307     1.0774056      1.214317
    ## 5      1.891680             0      1.723578     1.5698489      2.368597
    ##   0610010B08Rik;Gm14434;Gm14308
    ## 1                             0
    ## 2                             0
    ## 3                             0
    ## 4                             0
    ## 5                             0

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
