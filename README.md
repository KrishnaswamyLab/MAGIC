Markov Affinity-based Graph Imputation of Cells (MAGIC)
-------------------------------------------------------
van Dijk, David, et al. "MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data." BioRxiv (2017): 111591.

http://www.biorxiv.org/content/early/2017/02/25/111591

MAGIC has been implemented in Python3, Matlab, and R.

## Python3

#### Installation and dependencies for the Python version
1. The Python3 version of MAGIC can be installed using:

        $> git clone git://github.com/KrishnaswamyLab/magic.git
        $> cd magic
        $> sudo -H pip3 install .

2. MAGIC depends on a number of `python3` packages available on pypi and these dependencies are listed in `setup.py`
All the dependencies will be automatically installed using the above commands

3. After pulling updates to MAGIC from github, the package must be uninstalled and reinstalled:
		
		$> sudo -H pip3 uninstall magic
		$> sudo -H pip3 install .
		
#### Usage

##### Interactive command line
We have included two tutorial notebooks on MAGIC usage and results visualization for single cell RNA-seq data.

EMT data notebook: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/develop/python/tutorial_notebooks/Magic_single_cell_RNAseq_EMT_data.ipynb

Bone Marrow data notebook: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/develop/python/tutorial_notebooks/Magic_single_cell_RNAseq_bone_marrow_data.ipynb

##### GUI
A python GUI is now available for MAGIC. After following the installation steps listed below, the GUI can be invoked using

        $> magic_gui.py

##### Command line script
MAGIC can be run using the command line script `MAGIC.py` with the following parameters:

		$> MAGIC.py -h
		usage: MAGIC.py [-h] -d D -o O [-g G] [--gene-name-file GN]
                [--use-ensemble-ids] [--cell-axis CA] [--skip-rows SKIP_ROWS]
                [--skip-columns SKIP_COLUMNS] [-n] [-l L]
                [--mols-per-cell-min MOLS_PER_CELL_MIN]
                [--mols-per-cell-max MOLS_PER_CELL_MAX] [-p P]
                [--pca-non-random] [-t T] [-k K] [-ka KA] [-e E] [-r R]
                [--plot] [--t-max TM] [--n-genes NG]
                {csv,10x,10x_HDF5,mtx}

		run MAGIC

		positional arguments:
		  {csv,10x,10x_HDF5,mtx}
		                        what is the file type of your original data?

		optional arguments:
		  -h, --help            show this help message and exit

		data loading parameters:
		  -d D, --data-file D   File path of input data file.
		  -o O, --output-file O
		                        File path of where to save the MAGIC imputed data (in
		                        csv format).
		  -g G, --genome G      Genome must be specified when loading 10x_HDF5 data.
		  --gene-name-file GN   Gene name file must be specified when loading mtx
		                        data.
		  --use-ensemble-ids    Use ensemble IDs instead of gene names.
		  --cell-axis CA        When loading a csv, specify whether cells are on rows
		                        or columns (Default = 'rows').
		  --skip-rows SKIP_ROWS
		                        When loading a csv, number of rows to skip after the
		                        header row (Default = 0).
		  --skip-columns SKIP_COLUMNS
		                        When loading a csv, number of columns to skip after
		                        the header columns (Default = 0).

		normalization/filtering parameters:
		  -n, --no-normalize    Do not perform library size normalization on the data
		  -l L, --log-transform L
		                        Log-transform data with the specified pseudocount.
		  --mols-per-cell-min MOLS_PER_CELL_MIN
		                        Minimum molecules/cell to use in filtering.
		  --mols-per-cell-max MOLS_PER_CELL_MAX
		                        Maximum molecules/cell to use in filtering.

		MAGIC parameters:
		  -p P, --pca-components P
		                        Number of pca components to use when running MAGIC
		                        (Default = 20).
		  --pca-non-random      Do not used randomized solver in PCA computation.
		  -t T                  t parameter for running MAGIC. Default = None, in this
		                        case, the optimal t will be calculated .
		  -k K                  Number of nearest neighbors to use when running MAGIC
		                        (Default = 30).
		  -ka KA                knn-autotune parameter for running MAGIC (Default =
		                        10).
		  -e E, --epsilon E     Epsilon parameter for running MAGIC (Default = 1).
		  -r R, --rescale R     Percentile to rescale data to after running MAGIC
		                        (Default = 99).
		  --plot                Plot R2 plot generated in optimal t calculation
		                        (Default=False).
		  --t-max TM            Maximum t value used in optimal t calculation
		                        (Default=12).
		  --n-genes NG          Number of genes to use in optimal t calculation, a
		                        smaller number of genes speeds up the calculation
		                        (Default=500).
                        
## Matlab

#### Instructions for the Matlab version
1. run_magic.m -- MAGIC imputation function
2. test_magic.m -- Shows how to run MAGIC. Also included is a function for loading 10x format data (load_10x.m)

## R

#### Installation and dependencies for the R version
1. The R version of MAGIC can be installed using:

        $> library("devtools")
        $> install_github("KrishnaswamyLab/magic/R")

2. MAGIC depends on a number of `R` packages and these dependencies are listed in `DESCRIPTION`
All the dependencies will be automatically installed using the above commands

#### Usage

After installing the package, MAGIC can be run by loading the library and calling `run_magic()`:
	
	$> library(Rmagic)
	$> MAGIC_data <- run_magic(data, t=6, rescale_percent=0.99)
For a working example, see `R/tests/test_magic.R`. Please unzip the data provided in the `data` folder.
