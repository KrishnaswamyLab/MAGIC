# author: Scott Gigante <scott.gigante@yale.edu>
# (C) 2017 Krishnaswamy Lab GPLv2

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rc, animation
import numbers
import scprep
from scipy import sparse

from .magic import MAGIC
from .utils import in_ipynb


def _validate_gene(gene, data):
    if isinstance(gene, str):
        if not isinstance(data, pd.DataFrame):
            raise ValueError(
                "Non-integer gene names only valid with pd.DataFrame "
                "input. X is a {}, gene = {}".format(type(data).__name__, gene)
            )
        if gene not in data.columns:
            raise ValueError("gene {} not found".format(gene))
    elif gene is not None and not isinstance(gene, numbers.Integral):
        raise TypeError("Expected int or str. Got {}".format(type(gene).__name__))
    return gene


def animate_magic(
    data,
    gene_x,
    gene_y,
    gene_color=None,
    t_max=20,
    delay=2,
    operator=None,
    filename=None,
    ax=None,
    figsize=None,
    s=1,
    cmap="inferno",
    interval=200,
    dpi=100,
    ipython_html="jshtml",
    verbose=False,
    **kwargs
):
    """Animate a gene-gene relationship with increased diffusion

    Parameters
    ----------
    data: array-like
        Input data matrix
    gene_x : int or str
        Gene to put on the x axis
    gene_y : int or str
        Gene to put on the y axis
    gene_color : int or str, optional (default: None)
        Gene to color by. If None, no color vector is used
    t_max : int, optional (default: 20)
        maximum value of t to include in the animation
    delay : int, optional (default: 5)
        number of frames to dwell on the first frame before applying MAGIC
    operator : magic.MAGIC, optional (default: None)
        precomputed MAGIC operator. If None, one is created.
    filename : str, optional (default: None)
        If not None, saves a .gif or .mp4 with the output
    ax : `matplotlib.Axes` or None, optional (default: None)
        axis on which to plot. If None, an axis is created
    figsize : tuple, optional (default: None)
        Tuple of floats for creation of new `matplotlib` figure. Only used if
        `ax` is None.
    s : int, optional (default: 1)
        Point size
    cmap : str or callable, optional (default: 'inferno')
        Matplotlib colormap
    interval : float, optional (default: 30)
        Time in milliseconds between frames
    dpi : int, optional (default: 100)
        Dots per inch (image quality) in saved animation)
    ipython_html : {'html5', 'jshtml'}
        which html writer to use if using a Jupyter Notebook
    verbose : bool, optional (default: False)
        MAGIC operator verbosity
    *kwargs : arguments for MAGIC

    Returns
    -------
    A Matplotlib animation showing diffusion of an edge with increased t
    """
    if in_ipynb():
        # credit to
        # http://tiao.io/posts/notebooks/save-matplotlib-animations-as-gifs/
        rc("animation", html=ipython_html)

    if filename is not None:
        if filename.endswith(".gif"):
            writer = "imagemagick"
        elif filename.endswith(".mp4"):
            writer = "ffmpeg"
        else:
            raise ValueError(
                "filename must end in .gif or .mp4. Got {}".format(filename)
            )

    if operator is None:
        operator = MAGIC(verbose=verbose, **kwargs).fit(data)
    else:
        operator.set_params(verbose=verbose, **kwargs)
    gene_x = _validate_gene(gene_x, data)
    gene_y = _validate_gene(gene_y, data)
    gene_color = _validate_gene(gene_color, data)
    if gene_color is not None:
        genes = np.array([gene_x, gene_y, gene_color])
    else:
        genes = np.array([gene_x, gene_y])

    if isinstance(cmap, str):
        cmap = plt.cm.cmap_d[cmap]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        show = True
    else:
        fig = ax.get_figure()
        show = False

    data_magic = scprep.select.select_cols(data, idx=genes)
    data_magic = scprep.utils.toarray(data_magic)
    c = data_magic[gene_color] if gene_color is not None else None
    sc = ax.scatter(data_magic[gene_x], data_magic[gene_y], c=c, cmap=cmap)
    ax.set_title("t = 0")
    ax.set_xlabel(gene_x)
    ax.set_ylabel(gene_y)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if gene_color is not None:
        plt.colorbar(sc, label=gene_color, ticks=[])

    data_magic = [data]
    for t in range(t_max):
        operator.set_params(t=t + 1)
        data_magic.append(operator.transform(genes=genes))

    def init():
        return ax

    def animate(i):
        i = max(i - delay, 0)
        data_t = data_magic[i]
        data_t = data_t if isinstance(data, pd.DataFrame) else data_t.T
        sc.set_offsets(np.array([data_t[gene_x], data_t[gene_y]]).T)
        xlim = np.min(data_t[gene_x]), np.max(data_t[gene_x])
        xrange = xlim[1] - xlim[0]
        ax.set_xlim(xlim[0] - xrange / 10, xlim[1] + xrange / 10)
        ylim = np.min(data_t[gene_y]), np.max(data_t[gene_y])
        yrange = ylim[1] - ylim[0]
        ax.set_ylim(ylim[0] - yrange / 10, ylim[1] + yrange / 10)
        ax.set_title("t = {}".format(i))
        if gene_color is not None:
            color_t = data_t[gene_color]
            color_t -= np.min(color_t)
            color_t /= np.max(color_t)
            sc.set_facecolor(cmap(color_t))
        return ax

    ani = animation.FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=range(t_max + delay + 1),
        interval=interval,
        blit=False,
    )

    if filename is not None:
        ani.save(filename, writer=writer, dpi=dpi)

    if in_ipynb():
        # credit to https://stackoverflow.com/a/45573903/3996580
        plt.close()
    elif show:
        plt.tight_layout()
        fig.show()

    return ani
