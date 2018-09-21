# author: Scott Gigante <scott.gigante@yale.edu>
# (C) 2018 Krishnaswamy Lab GPLv2

from __future__ import print_function, division
import warnings
import scprep


def load_csv(filename, cell_axis='row', delimiter=',',
             gene_names=True, cell_names=True,
             sparse=False, **kwargs):
    """Load a csv file

    Parameters
    ----------
    filename : str
        The name of the csv file to be loaded
    cell_axis : {'row', 'column'}, optional (default: 'row')
        If your data has genes on the rows and cells on the columns, use
        cell_axis='column'
    delimiter : str, optional (default: ',')
        Use '\\t' for tab separated values (tsv)
    gene_names : `bool`, `str`, array-like, or `None` (default: True)
        If `True`, we assume gene names are in the first row/column. Otherwise
        expects a filename or an array containing a list of gene symbols or ids
    cell_names : `bool`, `str`, array-like, or `None` (default: True)
        If `True`, we assume cell names are in the first row/column. Otherwise
        expects a filename or an array containing a list of cell barcodes.
    sparse : bool, optional (default: False)
        If True, loads the data as a pd.SparseDataFrame. This uses less memory
        but more CPU.
    **kwargs : optional arguments for `pd.read_csv`.

    Returns
    -------
    data : pd.DataFrame
    """
    warnings.warn("magic.io is deprecated. Please use scprep.io instead. "
                  "Read more at http://scprep.readthedocs.io",
                  FutureWarning)
    return scprep.io.load_csv(filename=filename, cell_axis=cell_axis,
                              delimiter=delimiter,
                              gene_names=gene_names, cell_names=cell_names,
                              sparse=sparse, **kwargs)


def load_tsv(filename, cell_axis='row', delimiter='\t',
             gene_names=True, cell_names=True,
             sparse=False, **kwargs):
    """Load a tsv file

    Parameters
    ----------
    filename : str
        The name of the csv file to be loaded
    cell_axis : {'row', 'column'}, optional (default: 'row')
        If your data has genes on the rows and cells on the columns, use
        cell_axis='column'
    delimiter : str, optional (default: '\\t')
        Use ',' for comma separated values (csv)
    gene_names : `bool`, `str`, array-like, or `None` (default: True)
        If `True`, we assume gene names are in the first row/column. Otherwise
        expects a filename or an array containing a list of gene symbols or ids
    cell_names : `bool`, `str`, array-like, or `None` (default: True)
        If `True`, we assume cell names are in the first row/column. Otherwise
        expects a filename or an array containing a list of cell barcodes.
    sparse : bool, optional (default: False)
        If True, loads the data as a pd.SparseDataFrame. This uses less memory
        but more CPU.
    **kwargs : optional arguments for `pd.read_csv`.

    Returns
    -------
    data : pd.DataFrame
    """
    return load_csv(filename, cell_axis=cell_axis, delimiter=delimiter,
                    gene_names=gene_names, cell_names=cell_names,
                    sparse=sparse, **kwargs)


def load_fcs(filename, gene_names=True, cell_names=True,
             sparse=None,
             metadata_channels=['Time', 'Event_length', 'DNA1', 'DNA2',
                                'Cisplatin', 'beadDist', 'bead1']):
    """Load a fcs file

    Parameters
    ----------
    filename : str
        The name of the fcs file to be loaded
    gene_names : `bool`, `str`, array-like, or `None` (default: True)
        If `True`, we assume gene names are contained in the file. Otherwise
        expects a filename or an array containing a list of gene symbols or ids
    cell_names : `bool`, `str`, array-like, or `None` (default: True)
        If `True`, we assume cell names are contained in the file. Otherwise
        expects a filename or an array containing a list of cell barcodes.
    sparse : bool, optional (default: None)
        If True, loads the data as a pd.SparseDataFrame. This uses less memory
        but more CPU.
    metadata_channels : list-like, optional (default: ['Time', 'Event_length', 'DNA1', 'DNA2', 'Cisplatin', 'beadDist', 'bead1'])
        Channels to be excluded from the data

    Returns
    -------
    data : pd.DataFrame
    """
    warnings.warn("magic.io is deprecated. Please use scprep.io instead. "
                  "Read more at http://scprep.readthedocs.io",
                  FutureWarning)
    return scprep.io.load_fcs(filename=filename, gene_names=gene_names,
                              cell_names=cell_names,
                              sparse=sparse,
                              metadata_channels=metadata_channels)


def load_mtx(mtx_file, cell_axis='row',
             gene_names=None, cell_names=None, sparse=None):
    """Load a mtx file

    Parameters
    ----------
    filename : str
        The name of the mtx file to be loaded
    cell_axis : {'row', 'column'}, optional (default: 'row')
        If your data has genes on the rows and cells on the columns, use
        cell_axis='column'
    gene_names : `str`, array-like, or `None` (default: None)
        Expects a filename or an array containing a list of gene symbols or ids
    cell_names : `str`, array-like, or `None` (default: None)
        Expects a filename or an array containing a list of cell barcodes.
    sparse : bool, optional (default: None)
        If True, loads the data as a pd.SparseDataFrame. This uses less memory
        but more CPU.

    Returns
    -------
    data : pd.DataFrame
    """
    warnings.warn("magic.io is deprecated. Please use scprep.io instead. "
                  "Read more at http://scprep.readthedocs.io",
                  FutureWarning)
    return scprep.io.load_mtx(mtx_file=mtx_file, cell_axis=cell_axis,
                              gene_names=gene_names, cell_names=cell_names,
                              sparse=sparse)


def load_10X(data_dir, sparse=True, gene_labels='symbol',
             allow_duplicates=None):
    """Basic IO for 10X data produced from the 10X Cellranger pipeline.

    A default run of the `cellranger count` command will generate gene-barcode
    matrices for secondary analysis. For both "raw" and "filtered" output,
    directories are created containing three files:
    'matrix.mtx', 'barcodes.tsv', 'genes.tsv'.
    Running `phate.io.load_10X(data_dir)` will return a Pandas DataFrame with
    genes as columns and cells as rows. The returned DataFrame will be ready to
    use with PHATE.

    Parameters
    ----------
    data_dir: string
        path to input data directory
        expects 'matrix.mtx', 'genes.tsv', 'barcodes.tsv' to be present and
        will raise an error otherwise
    sparse: boolean
        If True, a sparse Pandas DataFrame is returned.
    gene_labels: string, {'id', 'symbol', 'both'} optional, default: 'symbol'
        Whether the columns of the dataframe should contain gene ids or gene
        symbols. If 'both', returns symbols followed by ids in parentheses.
    allow_duplicates : bool, optional (default: None)
        Whether or not to allow duplicate gene names. If None, duplicates are
        allowed for dense input but not for sparse input.

    Returns
    -------
    data: pandas.DataFrame shape = (n_cell, n_genes)
        imported data matrix
    """
    warnings.warn("magic.io is deprecated. Please use scprep.io instead. "
                  "Read more at http://scprep.readthedocs.io",
                  FutureWarning)
    return scprep.io.load_10X(data_dir=data_dir, sparse=sparse,
                              gene_labels=gene_labels,
                              allow_duplicates=allow_duplicates)


def load_10X_zip(filename, sparse=True, gene_labels='symbol',
                 allow_duplicates=None):
    """Basic IO for zipped 10X data produced from the 10X Cellranger pipeline.

    Runs `load_10X` after unzipping the data contained in `filename`

    Parameters
    ----------
    filename: string
        path to zipped input data directory
        expects 'matrix.mtx', 'genes.tsv', 'barcodes.tsv' to be present and
        will raise an error otherwise
    sparse: boolean
        If True, a sparse Pandas DataFrame is returned.
    gene_labels: string, {'id', 'symbol', 'both'} optional, default: 'symbol'
        Whether the columns of the dataframe should contain gene ids or gene
        symbols. If 'both', returns symbols followed by ids in parentheses.
    allow_duplicates : bool, optional (default: None)
        Whether or not to allow duplicate gene names. If None, duplicates are
        allowed for dense input but not for sparse input.

    Returns
    -------
    data: pandas.DataFrame shape = (n_cell, n_genes)
        imported data matrix
    """
    return scprep.io.load_10X_zip(filename=filename, sparse=sparse,
                                  gene_labels=gene_labels,
                                  allow_duplicates=allow_duplicates)


def load_10X_HDF5(filename, genome=None, sparse=True, gene_labels='symbol',
                  allow_duplicates=None):
    """Basic IO for HDF5 10X data produced from the 10X Cellranger pipeline.

    Equivalent to `load_10X` but for HDF5 format.

    Parameters
    ----------
    filename: string
        path to HDF5 input data
    genome : str or None, optional (default: None)
        Name of the genome to which CellRanger ran analysis. If None, selects
        the first available genome, and prints all available genomes if more
        than one is available.
    sparse: boolean
        If True, a sparse Pandas DataFrame is returned.
    gene_labels: string, {'id', 'symbol', 'both'} optional, default: 'symbol'
        Whether the columns of the dataframe should contain gene ids or gene
        symbols. If 'both', returns symbols followed by ids in parentheses.
    allow_duplicates : bool, optional (default: None)
        Whether or not to allow duplicate gene names. If None, duplicates are
        allowed for dense input but not for sparse input.

    Returns
    -------
    data: array-like, shape=[n_samples, n_features]
        If sparse, data will be a pd.SparseDataFrame. Otherwise, data will
        be a pd.DataFrame.
    """
    warnings.warn("magic.io is deprecated. Please use scprep.io instead. "
                  "Read more at http://scprep.readthedocs.io",
                  FutureWarning)
    return scprep.io.load_10X_HDF5(filename=filename, genome=genome,
                                   sparse=sparse,
                                   gene_labels=gene_labels,
                                   allow_duplicates=allow_duplicates)
