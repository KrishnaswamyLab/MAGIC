Tutorial
--------

To run MAGIC on your dataset, create a MAGIC operator and run `fit_transform`. Here we show an example with an artificial test dataset located in the MAGIC repository::

    import magic
    import matplotlib.pyplot as plt
    import pandas as pd
    X = pd.read_csv("MAGIC/data/test_data.csv")
    magic_operator = magic.MAGIC()
    X_magic = magic_operator.fit_transform(X, genes=['VIM', 'CDH1', 'ZEB1'])
    plt.scatter(X_magic['VIM'], X_magic['CDH1'], c=X_magic['ZEB1'], s=1, cmap='inferno')
    plt.show()
    magic.plot.animate_magic(X, gene_x='VIM', gene_y='CDH1', gene_color='ZEB1', operator=magic_operator)

A demo on MAGIC usage for single cell RNA-seq data can be found in this notebook_: `http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/master/python/tutorial_notebooks/emt_tutorial.ipynb`__

.. _notebook: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/master/python/tutorial_notebooks/emt_tutorial.ipynb

__ notebook_

A second tutorial analyzing myeloid and erythroid cells in mouse bone marrow is available here_: `http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/master/python/tutorial_notebooks/bonemarrow_tutorial.ipynb`__

.. _here: http://nbviewer.jupyter.org/github/KrishnaswamyLab/magic/blob/master/python/tutorial_notebooks/bonemarrow_tutorial.ipynb

__ here_
