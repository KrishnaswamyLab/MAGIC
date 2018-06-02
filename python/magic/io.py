# author: Daniel Burkhardt <daniel.burkhardt@yale.edu>
# (C) 2017 Krishnaswamy Lab GPLv2

from __future__ import print_function, division
import pandas as pd
import scipy.io as sio
import warnings
import numpy as np


def load_10X(data_dir, sparse=True, gene_labels='symbol'):
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
    data_dir : string
        path to input data directory
        expects 'matrix.mtx', 'genes.tsv', 'barcodes.csv' to be present and
        will raise and error otherwise
    sparse : boolean
        If True, a sparse Pandas DataFrame is returned.
    gene_labels : string, 'id' or 'symbol', optional, default: 'id'
        Whether the columns of the dataframe should contain gene ids or gene
        symbols

    Returns
    -------
    data : pandas.DataFrame shape=(n_cell, n_genes)
        imported data matrix
    """

    if gene_labels not in ['id', 'symbol']:
        raise ValueError("gene_labels not in ['id', 'symbol']")

    try:
        m = sio.mmread(data_dir + "/matrix.mtx")
        genes = pd.read_csv(data_dir + "/genes.tsv",
                            delimiter='\t', header=None)
        genes.columns = pd.Index(['id', 'symbol'])
        barcodes = pd.read_csv(data_dir + "/barcodes.tsv",
                               delimiter='\t', header=None)

    except (FileNotFoundError, OSError):
        raise FileNotFoundError(
            "'matrix.mtx', 'genes.tsv', and 'barcodes.tsv' must be present "
            "in data_dir")

    index = pd.Index(barcodes[0])
    columns = pd.Index(genes[gene_labels])
    if sparse and np.sum(columns.duplicated()) > 0:
        warnings.warn("Duplicate gene names detected! Forcing dense matrix. "
                      "Alternatively, try loading the matrix with "
                      "`gene_labels='id'`", RuntimeWarning)
        sparse = False

    if sparse:
        data = pd.SparseDataFrame(m.T, index=index,
                                  columns=columns,
                                  default_fill_value=0)
    else:
        data = pd.DataFrame(m.toarray().T, index=index, columns=columns)

    print("Imported data matrix with %s cells and %s genes." %
          (data.shape[0], data.shape[1]))
    return data

# TODO: clean these up.
    @classmethod
    def from_csv(cls, counts_csv_file, data_type='sc-seq', cell_axis=0, delimiter=',',
                 rows_after_header_to_skip=0, cols_after_header_to_skip=0, normalize=True):

        if not data_type in ['sc-seq', 'masscyt']:
            raise RuntimeError('data_type must be either sc-seq or masscyt')

        # Read in csv file
        df = pd.DataFrame.from_csv(counts_csv_file, sep=delimiter)

        df.drop(df.index[1:rows_after_header_to_skip + 1],
                axis=0, inplace=True)
        df.drop(df.columns[1:cols_after_header_to_skip + 1],
                axis=1, inplace=True)

        if cell_axis != 0:
            df = df.transpose()

        # Construct class object
        scdata = cls(df, data_type=data_type)

        # Normalize if specified
        if normalize == True:
            scdata = scdata.normalize_scseq_data()

        return scdata

    @classmethod
    def from_fcs(cls, fcs_file, cofactor=5,
                 metadata_channels=['Time', 'Event_length', 'DNA1', 'DNA2', 'Cisplatin', 'beadDist', 'bead1']):
        import fcsparser

        # Parse the fcs file
        text, data = fcsparser.parse(fcs_file)
        data = data.astype(np.float64)

        # Extract the S and N features (Indexing assumed to start from 1)
        # Assumes channel names are in S
        no_channels = text['$PAR']
        channel_names = [''] * no_channels
        for i in range(1, no_channels + 1):
            # S name
            try:
                channel_names[i - 1] = text['$P%dS' % i]
            except KeyError:
                channel_names[i - 1] = text['$P%dN' % i]
        data.columns = channel_names

        # Metadata and data
        metadata_channels = data.columns.intersection(metadata_channels)
        data_channels = data.columns.difference(metadata_channels)
        metadata = data[metadata_channels]
        data = data[data_channels]

        # Transform if necessary
        if cofactor is not None or cofactor > 0:
            data = np.arcsinh(np.divide(data, cofactor))

        # Create and return scdata object
        scdata = cls(data, 'masscyt', metadata)
        return scdata

    @classmethod
    def from_mtx(cls, mtx_file, gene_name_file, normalize=True):

        # Read in mtx file
        count_matrix = mmread(mtx_file)

        gene_names = np.loadtxt(gene_name_file, dtype=np.dtype('S'))
        gene_names = np.array([gene.decode('utf-8') for gene in gene_names])

        # remove todense
        df = pd.DataFrame(count_matrix.todense(), columns=gene_names)

        # Construct class object
        scdata = cls(df, data_type='sc-seq')

        # Normalize if specified
        if normalize == True:
            scdata = scdata.normalize_scseq_data()

        return scdata

    @classmethod
    def from_10x(cls, data_dir, use_ensemble_id=True, normalize=True):
        # loads 10x sparse format data
        # data_dir is dir that contains matrix.mtx, genes.tsv and barcodes.tsv
        # return_sparse=True -- returns data matrix in sparse format (default =
        # False)

        if data_dir == None or len(data_dir) == 0:
            data_dir = './'
        elif data_dir[len(data_dir) - 1] != '/':
            data_dir = data_dir + '/'

        filename_dataMatrix = os.path.expanduser(data_dir + 'matrix.mtx')
        filename_genes = os.path.expanduser(data_dir + 'genes.tsv')
        filename_cells = os.path.expanduser(data_dir + 'barcodes.tsv')

        # Read in gene expression matrix (sparse matrix)
        #Rows = genes, columns = cells
        print('LOADING')
        dataMatrix = mmread(filename_dataMatrix)

        # Read in row names (gene names / IDs)
        gene_names = np.loadtxt(
            filename_genes, delimiter='\t', dtype=bytes).astype(str)
        if use_ensemble_id == True:
            gene_names = [gene[0] for gene in gene_names]
        else:
            gene_names = [gene[1] for gene in gene_names]
        cell_names = np.loadtxt(
            filename_cells, delimiter='\t', dtype=bytes).astype(str)

        dataMatrix = pd.DataFrame(
            dataMatrix.todense(), columns=cell_names, index=gene_names)

        # combine duplicate genes
        if use_ensemble_id == False:
            dataMatrix = dataMatrix.groupby(dataMatrix.index).sum()

        dataMatrix = dataMatrix.transpose()

        # Remove empty cells
        print('Removing empty cells')
        cell_sums = dataMatrix.sum(axis=1)
        to_keep = np.where(cell_sums > 0)[0]
        dataMatrix = dataMatrix.ix[
            dataMatrix.index[to_keep], :].astype(np.float32)

        # Remove empty genes
        print('Removing empty genes')
        gene_sums = dataMatrix.sum(axis=0)
        to_keep = np.where(gene_sums > 0)[0]
        dataMatrix = dataMatrix.ix[:, to_keep].astype(np.float32)

        # Construct class object
        scdata = cls(dataMatrix, data_type='sc-seq')

        # Normalize if specified
        if normalize == True:
            scdata = scdata.normalize_scseq_data()

        return scdata

    @classmethod
    def from_10x_HDF5(cls, filename, genome, use_ensemble_id=True, normalize=True):
        import tables

        with tables.open_file(filename, 'r') as f:
            try:
                group = f.get_node(f.root, genome)
            except tables.NoSuchNodeError:
                print("That genome does not exist in this file.")
                return None
            if use_ensemble_id:
                gene_names = getattr(group, 'genes').read()
            else:
                gene_names = getattr(group, 'gene_names').read()
            barcodes = getattr(group, 'barcodes').read()
            data = getattr(group, 'data').read()
            indices = getattr(group, 'indices').read()
            indptr = getattr(group, 'indptr').read()
            shape = getattr(group, 'shape').read()
            matrix = csc_matrix((data, indices, indptr), shape=shape)

            dataMatrix = pd.DataFrame(matrix.todense(), columns=np.array([b.decode() for b in barcodes]),
                                      index=np.array([g.decode() for g in gene_names]))
            dataMatrix = dataMatrix.transpose()

            # Remove empty cells
            print('Removing empty cells')
            cell_sums = dataMatrix.sum(axis=1)
            to_keep = np.where(cell_sums > 0)[0]
            dataMatrix = dataMatrix.ix[
                dataMatrix.index[to_keep], :].astype(np.float32)

            # Remove empty genes
            print('Removing empty genes')
            gene_sums = dataMatrix.sum(axis=0)
            to_keep = np.where(gene_sums > 0)[0]
            dataMatrix = dataMatrix.ix[:, to_keep].astype(np.float32)

            # Construct class object
            scdata = cls(dataMatrix, data_type='sc-seq')

            # Normalize if specified
            if normalize == True:
                scdata = scdata.normalize_scseq_data()

            return scdata
