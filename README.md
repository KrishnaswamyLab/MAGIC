Markov Affinity-based Graph Imputation of Cells (MAGIC)
-------------------------------------------------------
van Dijk, David, et al. "MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data." BioRxiv (2017): 111591.

http://www.biorxiv.org/content/early/2017/02/25/111591

MAGIC has been implemented in Python3, Matlab, and R.

<p align="center">
<img src="https://github.com/KrishnaswamyLab/MAGIC/blob/master/magic.gif"/>
<br>
<i>Magic reveals the interaction between Vimentin (VIM), Cadherin-1 (CDH1), and Zinc finger E-box-binding homeobox 1 (ZEB1, encoded by colors).
</i>
</p>

## Python

### Installation

#### Installation with pip

To install with `pip`, run the following from a terminal:

        pip install --user git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python

#### Installation from GitHub

To clone the repository and install manually, run the following from a terminal:

        git clone git://github.com/KrishnaswamyLab/MAGIC.git
        cd MAGIC/python
        python setup.py install --user

### Usage

##### Interactive command line
We have included two tutorial notebooks on MAGIC usage and results visualization for single cell RNA-seq data.

EMT data notebook: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/develop/python/tutorial_notebooks/Magic_single_cell_RNAseq_EMT_data.ipynb

Bone Marrow data notebook: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/develop/python/tutorial_notebooks/Magic_single_cell_RNAseq_bone_marrow_data.ipynb

## Matlab

#### Instructions for the Matlab version
1. run_magic.m -- MAGIC imputation function
2. test_magic.m -- Shows how to run MAGIC. Also included is a function for loading 10x format data (load_10x.m)

## R

### Installation

To use MAGIC, you will need to install both the R and Python packages.

#### Installation with `devtools` and `pip`

You can install `Rmagic` with `devtools` by running the following in R:

        if (!require(devtools)) install.packages("devtools")
        library(devtools)
        install_github("KrishnaswamyLab/magic/R")

You then need to install MAGIC in Python with `pip` by running the following from a terminal:

		pip install --user git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python

If `python` or `pip` are not installed, you will need to install them. We recommend [Miniconda3](https://conda.io/miniconda.html) to install Python and `pip` together, or otherwise you can install `pip` from https://pip.pypa.io/en/stable/installing/.

#### Installation from GitHub

To clone the repository and install manually, run the following from a terminal:

        git clone git://github.com/KrishnaswamyLab/MAGIC.git
        cd MAGIC/python
        python setup.py install --user
        cd ../Rmagic
        R CMD INSTALL .

#### Usage

After installing the package, MAGIC can be run by loading the library and calling `run_magic()`:

		library(Rmagic)
		data(magic_testdata)
		MAGIC_data <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"), rescale_percent=99)
		library(ggplot2)
		ggplot(MAGIC_data) +
		  geom_point(aes(x=VIM, y=CDH1, color=ZEB1))

For a working example, see the Rmarkdown tutorials at https://github.com/KrishnaswamyLab/MAGIC/blob/R_magic2/R/inst/examples/bonemarrow_tutorial.md and https://github.com/KrishnaswamyLab/MAGIC/blob/R_magic2/R/inst/examples/EMT_tutorial.md or in `R/inst/examples`.
