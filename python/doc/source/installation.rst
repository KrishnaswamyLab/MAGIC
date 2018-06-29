Installation
============

Python installation
-------------------

Installation with `pip`
~~~~~~~~~~~~~~~~~~~~~~~

The Python version of MAGIC can be installed using::

        pip install --user git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python

Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

The Python version of MAGIC can be installed from GitHub by running the following from a terminal::

       git clone --recursive git://github.com/KrishnaswamyLab/MAGIC.git
       cd MAGIC/python
       python setup.py install --user

MATLAB installation
-------------------

1. The MATLAB version of MAGIC can be accessed using::

    git clone git://github.com/KrishnaswamyLab/MAGIC.git
    cd MAGIC/Matlab

2. Add the MAGIC/Matlab directory to your MATLAB path and run any of our `run` or `test` scripts to get a feel for MAGIC.

R installation
--------------

In order to use MAGIC in R, you must also install the Python package.

Installation with `devtools` and `pip`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install `Rmagic` from CRAN by running the following code in R::

        if (!require(devtools)) install.packages("devtools")
        library(devtools)
        install_github("KrishnaswamyLab/magic/R")

Install `magic` in Python by running the following code from a terminal::

        pip install --user git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python

If `python` or `pip` are not installed, you will need to install them. We recommend Miniconda3_ to install Python and `pip` together, or otherwise you can install `pip` from https://pip.pypa.io/en/stable/installing/.

.. _Miniconda3: https://conda.io/miniconda.html)

Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

The latest source version of MAGIC can be accessed by running the following in a terminal::

    git clone https://github.com/KrishnaswamyLab/MAGIC.git
    cd MAGIC/Rmagic
    R CMD INSTALL
    cd ../python
    python setup.py install --user

If the `Rmagic` folder is empty, you have may forgotten to use the `--recursive` option for `git clone`. You can rectify this by running the following in a terminal::

    cd MAGIC
    git submodule init
    git submodule update
    cd Rmagic
    R CMD INSTALL
    cd ../python
    python setup.py install --user
