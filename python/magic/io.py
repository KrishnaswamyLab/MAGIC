# author: Scott Gigante <scott.gigante@yale.edu>
# (C) 2018 Krishnaswamy Lab GPLv2

from __future__ import print_function, division
import warnings
import scprep


def load_csv(
    filename,
    cell_axis="row",
    delimiter=",",
    gene_names=True,
    cell_names=True,
    sparse=False,
    **kwargs
):
    """magic.io is deprecated. Please use scprep.io instead.
    Read more at http://scprep.readthedocs.io/
    """
    raise RuntimeError(
        "magic.io is deprecated. Please use scprep.io instead. "
        "Read more at http://scprep.readthedocs.io",
        FutureWarning,
    )


def load_tsv(
    filename,
    cell_axis="row",
    delimiter="\t",
    gene_names=True,
    cell_names=True,
    sparse=False,
    **kwargs
):
    """magic.io is deprecated. Please use scprep.io instead.
    Read more at http://scprep.readthedocs.io/
    """
    raise RuntimeError(
        "magic.io is deprecated. Please use scprep.io instead. "
        "Read more at http://scprep.readthedocs.io",
        FutureWarning,
    )


def load_fcs(
    filename,
    gene_names=True,
    cell_names=True,
    sparse=None,
    metadata_channels=[
        "Time",
        "Event_length",
        "DNA1",
        "DNA2",
        "Cisplatin",
        "beadDist",
        "bead1",
    ],
):
    """magic.io is deprecated. Please use scprep.io instead.
    Read more at http://scprep.readthedocs.io/
    """
    raise RuntimeError(
        "magic.io is deprecated. Please use scprep.io instead. "
        "Read more at http://scprep.readthedocs.io",
        FutureWarning,
    )


def load_mtx(mtx_file, cell_axis="row", gene_names=None, cell_names=None, sparse=None):
    """magic.io is deprecated. Please use scprep.io instead.
    Read more at http://scprep.readthedocs.io/
    """
    raise RuntimeError(
        "magic.io is deprecated. Please use scprep.io instead. "
        "Read more at http://scprep.readthedocs.io",
        FutureWarning,
    )


def load_10X(data_dir, sparse=True, gene_labels="symbol", allow_duplicates=None):
    """magic.io is deprecated. Please use scprep.io instead.
    Read more at http://scprep.readthedocs.io/
    """
    raise RuntimeError(
        "magic.io is deprecated. Please use scprep.io instead. "
        "Read more at http://scprep.readthedocs.io",
        FutureWarning,
    )


def load_10X_zip(filename, sparse=True, gene_labels="symbol", allow_duplicates=None):
    """magic.io is deprecated. Please use scprep.io instead.
    Read more at http://scprep.readthedocs.io/
    """
    raise RuntimeError(
        "magic.io is deprecated. Please use scprep.io instead. "
        "Read more at http://scprep.readthedocs.io",
        FutureWarning,
    )


def load_10X_HDF5(
    filename, genome=None, sparse=True, gene_labels="symbol", allow_duplicates=None
):
    """magic.io is deprecated. Please use scprep.io instead.
    Read more at http://scprep.readthedocs.io/
    """
    raise RuntimeError(
        "magic.io is deprecated. Please use scprep.io instead. "
        "Read more at http://scprep.readthedocs.io",
        FutureWarning,
    )
