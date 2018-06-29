=======================================================
MAGIC - Markov Affinity-based Graph Imputation of Cells
=======================================================

.. raw:: html

    <a href="https://travis-ci.com/KrishnaswamyLab/MAGIC"><img src="https://api.travis-ci.com/KrishnaswamyLab/magic.svg?branch=master" alt="Travis CI Build"></a>

.. raw:: html

    <a href="https://magic.readthedocs.io/"><img src="https://img.shields.io/readthedocs/magic.svg" alt="Read the Docs"></img></a>

.. raw:: html

    <a href="https://doi.org/10.1016/j.cell.2018.05.061"><img src="https://zenodo.org/badge/DOI/10.1016/j.cell.2018.05.061.svg" alt="Cell Publication DOI"></a>

.. raw:: html

    <a href="https://twitter.com/KrishnaswamyLab"><img src="https://img.shields.io/twitter/follow/KrishnaswamyLab.svg?style=social&label=Follow" alt="Twitter"></a>

.. raw:: html

    <a href="https://github.com/KrishnaswamyLab/MAGIC/"><img src="https://img.shields.io/github/stars/KrishnaswamyLab/MAGIC.svg?style=social&label=Stars" alt="GitHub stars"></a>

MAGIC is a tool that shares information across similar cells, via data diffusion, to denoise the cell count matrix and fill in missing transcripts. To see how MAGIC can be applied to single-cell RNA-seq, elucidating the epithelial-to-mesenchymal transition, read our `publication in Cell`_.

`David van Dijk, et al. Recovering Gene Interactions from Single-Cell Data Using Data Diffusion. 2018. Cell.`__

.. _`publication in Cell`: https://www.cell.com/cell/abstract/S0092-8674(18)30724-4

__ `publication in Cell`_

.. toctree::
    :maxdepth: 2

    installation
    tutorial
    api

Quick Start
===========

To run MAGIC on your dataset, create a MAGIC operator and run `fit_transform`. Here we show an example with a small, artificial dataset located in the MAGIC repository::

        import magic
        import pandas as pd
        import matplotlib.pyplot as plt
        X = pd.read_csv("MAGIC/data/test_data.csv")
        magic_operator = magic.MAGIC()
        X_magic = magic_operator.fit_transform(X, genes=['VIM', 'CDH1', 'ZEB1'])
        plt.scatter(X_magic['VIM'], X_magic['CDH1'], c=X_magic['ZEB1'], s=1, cmap='inferno')
        plt.show()
        magic.plot.animate_magic(X, gene_x='VIM', gene_y='CDH1', gene_color='ZEB1', operator=magic_operator)

.. autoclass:: magic.MAGIC
    :members:
    :noindex:
