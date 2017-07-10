Markov Affinity-based Graph Imputation of Cells (MAGIC)
-------------------------------------------------------

MAGIC has been implemented in Python3 and Matlab.

#### Installation and dependencies for the Python version
1. The Python3 version of MAGIC can be installed using:

        $> git clone git://github.com/pkathail/magic.git
        $> cd magic
        $> sudo -H pip3 install .

2. MAGIC depends on a number of `python3` packages available on pypi and these dependencies are listed in `setup.py`
All the dependencies will be automatically installed using the above commands

3. After pulling updates to MAGIC from github, the package must be uninstalled and reinstalled:
		
		$> sudo -H pip3 uninstall magic
		$> sudo -H pip3 install .
		
#### Usage

##### Interactive command line
A tutorial on MAGIC usage and results visualization for single cell RNA-seq data can be found in this notebook: http://nbviewer.jupyter.org/github/pkathail/magic/blob/develop/notebooks/Magic_single_cell_RNAseq.ipynb


##### GUI
A python GUI is now available for MAGIC. After following the installation steps listed below, the GUI can be invoked using

        $> magic_gui.py

##### Command line script
MAGIC can be run using the command line script `MAGIC.py` with the following parameters:

		$> MAGIC.py -h
		
		usage: MAGIC.py [-h] -d I -o O [-n] [-l L] [-g G] [-p N] [--pca-non-random]
                		[-t T] [-k K] [-ka KA] [-e E] [-r R]
                		{csv,10x,mtx}

		run MAGIC

		positional arguments:
		  {csv,10x,mtx}         what is the file type of your original data?

		optional arguments:
		  -h, --help            show this help message and exit

		required arguments:
		  -d D, --data-file D   File path of input data file.
		  -o O, --output-file O
		                        File path of where to save the MAGIC imputed data (in
		                        csv format).

		normalization parameters:
		  -n, --no-normalize    Do not perform library size normalization on the data
		  -l L, --log-transform L
		                        Log-transform data with the specified pseudocount.

		MAGIC parameters:
		  -g G, --gene-name-file G
		                        Gene name file must be specified when loading mtx
		                        data.
		  -p P, --pca-components P
		                        Number of pca components to use when running MAGIC
		                        (Default = 20).
		  --pca-non-random      Do not used randomized solver in PCA computation.
		  -t T					t parameter for running MAGIC (Default = 6).
		  -k K					Number of nearest neighbors to use when running MAGIC
                        		(Default = 30).
		  -ka KA				knn-autotune parameter for running MAGIC (Default =
                        		10).
		  -e E, --epsilon E		Epsilon parameter for running MAGIC (Default = 1).
		  -r R, --rescale R		Percentile to rescale data to after running MAGIC
                        		(Default = 99).

#### Instructions for the Matlab version
1. run_magic.m -- MAGIC imputation function
2. test_magic.m -- Shows how to run MAGIC. Also included is a function for loading 10x format data (load_10x.m)
