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

#### Usage

##### Command line
A tutorial on MAGIC usage and results visualization for single cell RNA-seq data can be found in this notebook: http://nbviewer.jupyter.org/github/pkathail/magic/blob/develop/notebooks/Magic_single_cell_RNAseq.ipynb


##### GUI
A python GUI is now available for MAGIC. After following the installation steps listed below, the GUI can be invoked using

        $> magic_gui.py

#### Installation and dependencies for the Matlab version
1. run_magic.m uses Mauro Maggioni's Diffusion Geometry code. Download from here: http://www.math.jhu.edu/~mauro/Code/DiffusionGeometry_01.zip or use included DiffusionGeometry_01.zip
2. run_magic2.m is an implementation that does not use the Diffusion Geometry code.
3. test_magic.m shows how to run MAGIC. Also included is a function for loading 10x format data (load_10x.m)
