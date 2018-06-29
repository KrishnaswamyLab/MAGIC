Markov Affinity-based Graph Imputation of Cells (MAGIC)
-------------------------------------------------------

Markov Affinity-based Graph Imputation of Cells (MAGIC) is an algorithm for denoising and transcript recover of single cells applied to single-cell RNA sequencing data, as described in Van Dijk D *et al.* (2018), *Recovering Gene Interactions from Single-Cell Data Using Data Diffusion*, Cell <https://www.cell.com/cell/abstract/S0092-8674(18)30724-4>.

For R and MATLAB implementations of MAGIC, see <https://github.com/KrishnaswamyLab/MAGIC>.

<p align="center">
<img src="https://github.com/KrishnaswamyLab/MAGIC/blob/master/magic.gif"/>
<br>
<i>Magic reveals the interaction between Vimentin (VIM), Cadherin-1 (CDH1), and Zinc finger E-box-binding homeobox 1 (ZEB1, encoded by colors).
</i>
</p>

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

##### Example data

The following code runs MAGIC on test data located in the MAGIC repository.

		import magic
		import phate
		import pandas as pd
		X = pd.read_csv("MAGIC/data/test_data.csv")
		magic_operator = magic.MAGIC()
		X_magic = magic_operator.fit_transform(X, genes=['VIM', 'CDH1', 'ZEB1'])
		phate.plot.scatter2d(X_magic[['VIM', 'CDH1']], c=X_magic['ZEB1'])

##### Interactive command line
We have included two tutorial notebooks on MAGIC usage and results visualization for single cell RNA-seq data.

EMT data notebook: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/develop/python/tutorial_notebooks/Magic_single_cell_RNAseq_EMT_data.ipynb

Bone Marrow data notebook: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/develop/python/tutorial_notebooks/Magic_single_cell_RNAseq_bone_marrow_data.ipynb
