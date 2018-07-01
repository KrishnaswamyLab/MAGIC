# author: Scott Gigante <scott.gigante@yale.edu>
# (C) 2018 Krishnaswamy Lab GPLv2

from __future__ import print_function, division
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import warnings
import numpy as np
import os
import zipfile
import tempfile
import shutil
try:
    import fcsparser
except ImportError:
    pass
try:
    import tables
except ImportError:
    pass

try:
    FileNotFoundError
except NameError:
    # py2 compatibility
    FileNotFoundError = OSError


def with_fcsparser(fun):
    def wrapped_fun(*args, **kwargs):
        try:
            fcsparser
        except NameError:
            raise ImportError(
                "fcsparser not found. "
                "Please install it with e.g. `pip install --user fcsparser`")
        return fun(*args, **kwargs)
    return wrapped_fun


def with_tables(fun):
    def wrapped_fun(*args, **kwargs):
        try:
            tables
        except NameError:
            raise ImportError(
                "tables not found. "
                "Please install it with e.g. `pip install --user tables`")
        return fun(*args, **kwargs)
    return wrapped_fun


def _parse_header(header, n_expected, header_type="gene_names"):
    """
    Parameters
    ----------
    header : `str` filename, array-like or `None`

    n_expected : `int`
        Expected header length

    header_type : argument name for error printing
    """
    if header is None or header is False:
        return None
    elif isinstance(header, str):
        # treat as a file
        if header.endswith("tsv"):
            delimiter = "\t"
        else:
            delimiter = ","
        columns = pd.read_csv(header, delimiter=delimiter,
                              header=None).values.reshape(-1)
        if not len(columns) == n_expected:
            raise ValueError("Expected {} entries in {}. Got {}".format(
                n_expected, header, len(columns)))
    else:
        # treat as list
        columns = header
        if not len(columns) == n_expected:
            raise ValueError("Expected {} entries in {}. Got {}".format(
                n_expected, header_type, len(columns)))
    return columns


def _parse_gene_names(header, data):
    return _parse_header(header, data.shape[1],
                         header_type="gene_names")


def _parse_cell_names(header, data):
    return _parse_header(header, data.shape[0],
                         header_type="cell_names")


def _matrix_to_data_frame(data, gene_names=None, cell_names=None, sparse=None):
    """Return the optimal data type given data, gene names and cell names.

    Parameters
    ----------

    data : array-like

    gene_names : `str`, array-like or `None` (default: None)
        Either a filename or an array containing a list of gene symbols or ids.

    cell_names : `str`, array-like or `None` (default: None)
        Either a filename or an array containing a list of cell barcodes.

    sparse : `bool` or `None` (default: None)
        If not `None`, overrides default sparsity of the data.
    """
    if gene_names is None and cell_names is None and \
            not isinstance(data, pd.DataFrame):
        # just a matrix
        if sparse is not None:
            if sparse:
                if not sp.issparse(data):
                    # return scipy.sparse.csr_matrix
                    data = sp.csr_matrix(data)
            elif sp.issparse(data) and not sparse:
                # return numpy.ndarray
                data = data.toarray()
        else:
            # return data as is
            pass
        return data
    else:
        gene_names = _parse_gene_names(gene_names, data)
        cell_names = _parse_cell_names(cell_names, data)
        # dataframe with index and/or columns
        if sparse is None:
            # let the input data decide
            sparse = isinstance(data, pd.SparseDataFrame) or sp.issparse(data)
        if sparse and gene_names is not None and \
                len(np.unique(gene_names)) < len(gene_names):
            warnings.warn(
                "Duplicate gene names detected! Forcing dense matrix",
                RuntimeWarning)
            sparse = False
        if sparse:
            # return pandas.SparseDataFrame
            if isinstance(data, pd.DataFrame):
                if gene_names is not None:
                    data.columns = gene_names
                if cell_names is not None:
                    data.index = cell_names
                if not isinstance(data, pd.SparseDataFrame):
                    data = data.to_sparse(fill_value=0.0)
            else:
                data = pd.SparseDataFrame(data, default_fill_value=0.0,
                                          index=cell_names, columns=gene_names)
        else:
            # return pandas.DataFrame
            if isinstance(data, pd.DataFrame):
                if gene_names is not None:
                    data.columns = gene_names
                if cell_names is not None:
                    data.index = cell_names
                if isinstance(data, pd.SparseDataFrame):
                    data = data.to_dense()
            else:
                if sp.issparse(data):
                    data = data.toarray()
                data = pd.DataFrame(data, index=cell_names, columns=gene_names)
        return data


def _read_csv_sparse(filename, chunksize=1000000, fill_value=0.0, **kwargs):
    chunks = pd.read_csv(filename, chunksize=chunksize, **kwargs)
    data = pd.concat(chunk.to_sparse(fill_value=fill_value)
                     for chunk in chunks)
    return data


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
    if cell_axis not in ['row', 'column', 'col']:
        raise ValueError(
            "cell_axis {} not recognized. Expected 'row' or 'column'".format(
                cell_axis))

    if 'index_col' in kwargs:
        # override
        index_col = kwargs['index_col']
        cell_names = None
        del kwargs['index_col']
    elif cell_names is True:
        index_col = 0
        cell_names = None
    else:
        index_col = None

    if 'header' in kwargs:
        # override
        header = kwargs['header']
        del kwargs['header']
        gene_names = None
    elif gene_names is True:
        header = 0
        gene_names = None
    else:
        header = None

    # Read in csv file
    if sparse:
        read_fun = _read_csv_sparse
    else:
        read_fun = pd.read_csv
    data = read_fun(filename, delimiter=delimiter,
                    header=header, index_col=index_col,
                    **kwargs)

    if cell_axis in ['column', 'col']:
        data = data.T

    data = _matrix_to_data_frame(
        data, gene_names=gene_names,
        cell_names=cell_names, sparse=sparse)
    return data


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


@with_fcsparser
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
    if cell_names is True:
        cell_names = None
    if gene_names is True:
        gene_names = None
    # Parse the fcs file
    meta, data = fcsparser.parse(filename)
    metadata_channels = data.columns.intersection(metadata_channels)
    data_channels = data.columns.difference(metadata_channels)
    metadata = data[metadata_channels]
    data = data[data_channels]
    data = _matrix_to_data_frame(data, gene_names=gene_names,
                                 cell_names=cell_names, sparse=sparse)
    return metadata, data


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
    if cell_axis not in ['row', 'column', 'col']:
        raise ValueError(
            "cell_axis {} not recognized. Expected 'row' or 'column'".format(
                cell_axis))
    # Read in mtx file
    data = sio.mmread(mtx_file)
    if cell_axis in ['column', 'col']:
        data = data.T
    data = _matrix_to_data_frame(
        data, gene_names=gene_names,
        cell_names=cell_names, sparse=sparse)
    return data


def _combine_gene_id(symbols, ids):
    """Creates gene labels of the form SYMBOL (ID)

    Parameters
    ----------

    genes: pandas.DataFrame with columns['symbol', 'id']

    Returns
    -------

    pandas.Index with combined gene symbols and ids
    """
    columns = np.core.defchararray.add(
        np.array(symbols, dtype=str), ' (')
    columns = np.core.defchararray.add(
        columns, np.array(ids, dtype=str))
    columns = np.core.defchararray.add(columns, ')')
    return columns


def _parse_10x_genes(symbols, ids, gene_labels='symbol',
                     allow_duplicates=True):
    if gene_labels not in ['symbol', 'id', 'both']:
        raise ValueError("gene_labels='{}' not recognized. Choose from "
                         "['symbol', 'id', 'both']")
    if gene_labels == 'both':
        columns = _combine_gene_id(symbols, ids)
    if gene_labels == 'symbol':
        columns = symbols
        if not allow_duplicates and len(np.unique(columns)) < len(columns):
            warnings.warn(
                "Duplicate gene names detected! Forcing `gene_labels='id'`. "
                "Alternatively, try `gene_labels='both'`, "
                "`allow_duplicates=True`, or load the matrix"
                " with `sparse=False`", RuntimeWarning)
            gene_labels = 'id'
    if gene_labels == 'id':
        columns = ids
    return columns


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

    if gene_labels not in ['id', 'symbol', 'both']:
        raise ValueError("gene_labels not in ['id', 'symbol', 'both']")

    if not os.path.isdir(data_dir):
        raise FileNotFoundError(
            "{} is not a directory".format(data_dir))

    try:
        m = sio.mmread(os.path.join(data_dir, "matrix.mtx"))
        genes = pd.read_csv(os.path.join(data_dir, "genes.tsv"),
                            delimiter='\t', header=None)
        genes.columns = ['id', 'symbol']
        barcodes = pd.read_csv(os.path.join(data_dir, "barcodes.tsv"),
                               delimiter='\t', header=None)

    except (FileNotFoundError, OSError):
        raise FileNotFoundError(
            "'matrix.mtx', 'genes.tsv', and 'barcodes.tsv' must be present "
            "in {}".format(data_dir))

    cell_names = barcodes[0]
    if allow_duplicates is None:
        allow_duplicates = not sparse
    gene_names = _parse_10x_genes(genes['symbol'], genes['id'],
                                  gene_labels=gene_labels,
                                  allow_duplicates=allow_duplicates)

    data = _matrix_to_data_frame(m.T, cell_names=cell_names,
                                 gene_names=gene_names,
                                 sparse=sparse)
    return data


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
    tmpdir = tempfile.mkdtemp()
    with zipfile.ZipFile(filename) as handle:
        files = handle.namelist()
        if len(files) != 4:
            valid = False
        else:
            dirname = files[0].strip("/")
            subdir_files = [f.split("/")[-1] for f in files]
            if "barcodes.tsv" not in subdir_files:
                valid = False
            elif "genes.tsv" not in subdir_files:
                valid = False
            elif "matrix.mtx" not in subdir_files:
                valid = False
            else:
                valid = True
        if not valid:
            raise ValueError(
                "Expected a single zipped folder containing 'matrix.mtx', "
                "'genes.tsv', and 'barcodes.tsv'. Got {}".format(files))
        handle.extractall(path=tmpdir)
    data = load_10X(os.path.join(tmpdir, dirname))
    shutil.rmtree(tmpdir)
    return data


@with_tables
def load_10x_HDF5(filename, genome, sparse=True, gene_labels='symbol',
                  allow_duplicates=None):
    """Basic IO for HDF5 10X data produced from the 10X Cellranger pipeline.

    Equivalent to `load_10X` but for HDF5 format.

    Parameters
    ----------
    filename: string
        path to HDF5 input data
    genome : str
        Name of the genome to which CellRanger ran analysis
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
    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            raise ValueError(
                "Genome {} not found in {}.".format(genome, filename))
            # TODO: print available genomes.
        if allow_duplicates is None:
            allow_duplicates = not sparse
        gene_names = _parse_10x_genes(
            symbols=[g.decode() for g in getattr(group, 'gene_names').read()],
            ids=[g.decode() for g in getattr(group, 'gene').read()],
            gene_labels=gene_labels, allow_duplicates=allow_duplicates)
        cell_names = [b.decode() for b in getattr(group, 'barcodes').read()]
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        data = sp.csr_matrix((data, indices, indptr), shape=shape)
        data = _matrix_to_data_frame(data.T,
                                     gene_names=gene_names,
                                     cell_names=cell_names,
                                     sparse=sparse)
        return data
