===========================================================================
MAGIC - Potential of Heat-diffusion for Affinity-based Trajectory Embedding
===========================================================================

.. raw:: html

    <a href="https://pypi.org/project/magic/"><img src="https://img.shields.io/pypi/v/magic.svg" alt="Latest PyPi version"></a>

.. raw:: html

    <a href="https://cran.r-project.org/package=magicR"><img src="https://img.shields.io/cran/v/magicR.svg" alt="Latest CRAN version"></a>

.. raw:: html

    <a href="https://travis-ci.com/KrishnaswamyLab/MAGIC"><img src="https://api.travis-ci.com/KrishnaswamyLab/magic.svg?branch=master" alt="Travis CI Build"></a>

.. raw:: html

    <a href="https://magic.readthedocs.io/"><img src="https://img.shields.io/readthedocs/magic.svg" alt="Read the Docs"></img></a>

.. raw:: html

    <a href="https://www.biorxiv.org/content/early/2017/12/01/120378"><img src="https://zenodo.org/badge/DOI/10.1101/120378.svg" alt="bioRxiv Preprint"></a>

.. raw:: html

    <a href="https://twitter.com/KrishnaswamyLab"><img src="https://img.shields.io/twitter/follow/KrishnaswamyLab.svg?style=social&label=Follow" alt="Twitter"></a>

.. raw:: html

    <a href="https://github.com/KrishnaswamyLab/MAGIC/"><img src="https://img.shields.io/github/stars/KrishnaswamyLab/MAGIC.svg?style=social&label=Stars" alt="GitHub stars"></a>

MAGIC is a tool for visualizing high dimensional single-cell data with natural progressions or trajectories. MAGIC uses a novel conceptual framework for learning and visualizing the manifold inherent to biological systems in which smooth transitions mark the progressions of cells from one state to another. To see how MAGIC can be applied to single-cell RNA-seq datasets from hematopoietic stem cells, human embryonic stem cells, and bone marrow samples, check out our `preprint on BioRxiv`_.

`Kevin R. Moon, David van Dijk, Zheng Wang, et al. MAGIC: A Dimensionality Reduction Method for Visualizing Trajectory Structures in High-Dimensional Biological Data. 2017. BioRxiv.`__

.. _`preprint on BioRxiv`: https://www.biorxiv.org/content/early/2017/03/24/120378

__ `preprint on BioRxiv`_

.. toctree::
    :maxdepth: 2

    installation
    tutorial
    api

Quick Start
===========

To run MAGIC on your dataset, create a MAGIC operator and run `fit_transform`. Here we show an example with an artificial tree::

        import magic
        tree_data, tree_clusters = magic.tree.gen_dla()
        magic_operator = magic.MAGIC(k=15, t=100)
        tree_magic = magic_operator.fit_transform(tree_data)
        magic.plot.scatter2d(magic_operator, c=tree_clusters) # or magic.plot.scatter2d(tree_magic, c=tree_clusters)
        magic.plot.rotate_scatter3d(magic_operator, c=tree_clusters)

.. autoclass:: magic.MAGIC
    :members:
    :noindex:
