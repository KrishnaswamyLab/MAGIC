import re
import os
import random
import pickle
import warnings
import shlex
import shutil
from copy import deepcopy
from collections import defaultdict, Counter
from subprocess import call, Popen, PIPE
import glob

import numpy as np
import pandas as pd

import matplotlib
# try:
#     os.environ['DISPLAY']
# except KeyError:
#     matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # catch experimental ipython widget warning
    import seaborn as sns

# from tsne import bh_sne
from sklearn.manifold import TSNE
from sklearn.manifold.t_sne import _joint_probabilities, _joint_probabilities_nn
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from scipy.spatial.distance import squareform
from scipy.sparse import csr_matrix, find, vstack, hstack, issparse
from scipy.sparse.linalg import eigs
from numpy.linalg import norm
from scipy.stats import gaussian_kde
from scipy.io import mmread
from numpy.core.umath_tests import inner1d

import fcsparser

import magic

# set plotting defaults
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # catch experimental ipython widget warning
    sns.set(context="paper", style='ticks', font_scale=1.5, font='Bitstream Vera Sans')

matplotlib.rcParams['image.cmap'] = 'viridis'
size = 8


def qualitative_colors(n):
    """ Generalte list of colors
    :param n: Number of colors
    """
    return sns.color_palette('Set1', n)


def get_fig(fig=None, ax=None, figsize=[6.5, 6.5]):
    """fills in any missing axis or figure with the currently active one
    :param ax: matplotlib Axis object
    :param fig: matplotlib Figure object
    """
    if not fig:
        fig = plt.figure(figsize=figsize)
    if not ax:
        ax = plt.gca()
    return fig, ax

def density_2d(x, y):
    """return x and y and their density z, sorted by their density (smallest to largest)

    :param x:
    :param y:
    :return:
    """
    xy = np.vstack([np.ravel(x), np.ravel(y)])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    return np.ravel(x)[i], np.ravel(y)[i], np.arcsinh(z[i])

        
class SCData:

    def __init__(self, data, data_type='sc-seq', metadata=None):
        """
        Container class for single cell data
        :param data:  DataFrame of cells X genes representing expression
        :param data_type: Type of the data: Can be either 'sc-seq' or 'masscyt'
        :param metadata: None or DataFrame representing metadata about the cells
        """
        if not (isinstance(data, pd.DataFrame)):
            raise TypeError('data must be of type or DataFrame')
        if not data_type in ['sc-seq', 'masscyt']:
            raise RuntimeError('data_type must be either sc-seq or masscyt')
        if metadata is None:
            metadata = pd.DataFrame(index=data.index, dtype='O')
        cols = [np.array(['data' for i in range(data.shape[1])]), np.array(data.columns.values)]
        self._data = pd.DataFrame(data.values, index=data.index, columns=cols)
        self._metadata = metadata
        self._data_type = data_type
        self._normalized = False
        self._pca = None
        self._tsne = None
        self._diffusion_eigenvectors = None
        self._diffusion_eigenvalues = None
        self._diffusion_map_correlations = None
        self._magic = None
        self._normalized = False

        # Library size
        self._library_sizes = None

    def save(self, fout: str):  # -> None:
        """
        :param fout: str, name of archive to store pickled SCData data in. Should end
          in '.p'.
        :return: None
        """
        with open(fout, 'wb') as f:
            pickle.dump(vars(self), f)

    @classmethod
    def load(cls, fin):
        """

        :param fin: str, name of pickled archive containing SCData data
        :return: SCData
        """
        with open(fin, 'rb') as f:
            data = pickle.load(f)
        scdata = cls(data['_data'], data['_metadata'])
        del data['_data']
        del data['_metadata']
        for k, v in data.items():
            setattr(scdata, k[1:], v)
        return scdata

    def __repr__(self):
        c, g = self.data.shape
        _repr = ('SCData: {c} cells x {g} genes\n'.format(g=g, c=c))
        for k, v in sorted(vars(self).items()):
            if not (k == '_data'):
                _repr += '\n{}={}'.format(k[1:], 'None' if v is None else 'True')
        return _repr

    @property
    def data_type(self):
        return self._data_type

    @property
    def data(self):
        return self._data['data']

    @data.setter
    def data(self, item):
        if not (isinstance(item, pd.DataFrame)):
            raise TypeError('SCData.data must be of type DataFrame')
        cols = [np.array(['data' for i in range(item.shape[1])]), np.array(item.columns.values)]
        self._data = pd.DataFrame(item.values, index=item.index, columns=cols)

    @property 
    def extended_data(self):
        return self._data

    @extended_data.setter
    def extended_data(self, item):
        if not (isinstance(item, pd.DataFrame)):
            raise TypeError('SCData.extended_data must be of type DataFrame')
        self._data = item

    @property
    def metadata(self):
        return self._metadata

    @metadata.setter
    def metadata(self, item):
        if not isinstance(item, pd.DataFrame):
            raise TypeError('SCData.metadata must be of type DataFrame')
        self._metadata = item

    @property
    def pca(self):
        if self._pca != None:
            return self._data['PCA']
        return self._pca

    @pca.setter
    def pca(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.pca must be a dictionary of pd.DataFrame object')
        if item is None:
            del self._data['PCA']
            self._pca = None
        else:
            item.columns = [np.array(['PCA' for i in range(item.shape[1])]), item.columns.values]
            self._data = pd.concat([self._data, item], axis=1)
            self._pca = True


    @property
    def tsne(self):
        if self._tsne != None:
            return self._data['tSNE']
        return self._tsne

    @tsne.setter
    def tsne(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.tsne must be a pd.DataFrame object')
        if item is None:
            del self._data['tSNE']
            self._tsne = None
        else:
            item.columns = [np.array(['tSNE' for i in range(item.shape[1])]), item.columns.values]
            self._data = pd.concat([self._data, item], axis=1)
            self._tsne = True

    @property
    def diffusion_eigenvectors(self):
        if self._diffusion_eigenvectors != None:
            return self._data['diffusion_eigenvectors']
        return self._diffusion_eigenvectors

    @diffusion_eigenvectors.setter
    def diffusion_eigenvectors(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.diffusion_eigenvectors must be a pd.DataFrame object')
        if item is None:
            del self._data['diffusion_eigenvectors']
            self._diffusion_eigenvectors = None
        else:
            item.columns = [np.array(['diffusion_eigenvectors' for i in range(item.shape[1])]), item.columns.values]
            self._data = pd.concat([self._data, item], axis=1)
            self._diffusion_eigenvectors = True

    @property
    def diffusion_eigenvalues(self):
        return self._diffusion_eigenvalues

    @diffusion_eigenvalues.setter
    def diffusion_eigenvalues(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.diffusion_eigenvalues must be a pd.DataFrame object')
        self._diffusion_eigenvalues = item

    @property
    def diffusion_map_correlations(self):
        return self._diffusion_map_correlations

    @diffusion_map_correlations.setter
    def diffusion_map_correlations(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.diffusion_map_correlations must be a pd.DataFrame'
                            'object')
        self._diffusion_map_correlations = item

    @property
    def magic(self):
        return self._magic

    @magic.setter
    def magic(self, item):
        if not (isinstance(item, magic.mg.SCData) or item is None):
            raise TypeError('self.magic must be a SCData object')
        self._magic = item

    @property
    def library_sizes(self):
        return self._library_sizes

    @library_sizes.setter
    def library_sizes(self, item):
        if not (isinstance(item, pd.Series) or item is None):
            raise TypeError('self.library_sizes must be a pd.Series object')


    @classmethod
    def from_csv(cls, counts_csv_file, data_type, cell_axis=0, delimiter=',', 
                 rows_after_header_to_skip=0, cols_after_header_to_skip=0, normalize=True):

        if not data_type in ['sc-seq', 'masscyt']:
            raise RuntimeError('data_type must be either sc-seq or masscyt')

        # Read in csv file
        df = pd.DataFrame.from_csv( counts_csv_file, sep=delimiter)

        df.drop(df.index[1:rows_after_header_to_skip+1], axis=0, inplace=True)
        df.drop(df.columns[1:cols_after_header_to_skip+1],axis=1, inplace=True)
        
        if cell_axis != 0:
            df = df.transpose()
            
        # Construct class object
        scdata = cls( df, data_type=data_type )

        # Normalize if specified
        if normalize==True:
            scdata = scdata.normalize_scseq_data( )

        return scdata


    @classmethod
    def from_fcs(cls, fcs_file, cofactor=5, 
        metadata_channels=['Time', 'Event_length', 'DNA1', 'DNA2', 'Cisplatin', 'beadDist', 'bead1']):

        # Parse the fcs file
        text, data = fcsparser.parse( fcs_file )
        data = data.astype(np.float64)

        # Extract the S and N features (Indexing assumed to start from 1)
        # Assumes channel names are in S
        no_channels = text['$PAR']
        channel_names = [''] * no_channels
        for i in range(1, no_channels+1):
            # S name
            try:
                channel_names[i - 1] = text['$P%dS' % i]
            except KeyError:
                channel_names[i - 1] = text['$P%dN' % i]
        data.columns = channel_names
        
        # Metadata and data 
        metadata_channels = data.columns.intersection(metadata_channels)
        data_channels = data.columns.difference( metadata_channels )
        metadata = data[metadata_channels]
        data = data[data_channels]

        # Transform if necessary
        if cofactor is not None or cofactor > 0:
            data = np.arcsinh(np.divide( data, cofactor ))

        # Create and return scdata object
        scdata = cls(data, 'masscyt', metadata)
        return scdata


    @classmethod
    def from_mtx(cls, mtx_file, gene_name_file):

        #Read in mtx file
        count_matrix = mmread(mtx_file)

        gene_names = np.loadtxt(gene_name_file, dtype=np.dtype('S'))
        gene_names = np.array([gene.decode('utf-8') for gene in gene_names])

        ### remove todense
        df = pd.DataFrame(count_matrix.todense(), columns=gene_names)

        # Construct class object
        scdata = cls( df, data_type='sc-seq' )

        return scdata

    @classmethod
    def from_10x(cls, data_dir, use_ensemble_id=True):
        #loads 10x sparse format data
        #data_dir is dir that contains matrix.mtx, genes.tsv and barcodes.tsv
        #return_sparse=True -- returns data matrix in sparse format (default = False)

        if data_dir==None or len(data_dir)==0:
            data_dir = './'
        elif data_dir[len(data_dir)-1] != '/':
            data_dir = data_dir + '/'

        filename_dataMatrix = os.path.expanduser(data_dir + 'matrix.mtx')
        filename_genes = os.path.expanduser(data_dir + 'genes.tsv')
        filename_cells = os.path.expanduser(data_dir + 'barcodes.tsv')


        #Read in gene expression matrix (sparse matrix)
        #Rows = genes, columns = cells
        print('LOADING')
        dataMatrix = mmread(filename_dataMatrix);

        #Read in row names (gene names / IDs)
        gene_names = np.loadtxt(filename_genes, delimiter='\t', dtype=bytes).astype(str)
        if use_ensemble_id==True:
            gene_names = [gene[0] for gene in gene_names]
        else:
            gene_names = [gene[1] for gene in gene_names]
        cell_names = np.loadtxt(filename_cells, delimiter='\t', dtype=bytes).astype(str)
   
        dataMatrix = pd.DataFrame(dataMatrix.todense(), columns=cell_names, index=gene_names)

        #combine duplicate genes
        if use_ensemble_id == False:
            dataMatrix = dataMatrix.groupby(dataMatrix.index).sum()

        dataMatrix = dataMatrix.transpose()

        #Remove empty cells
        print('Removing empty cells')
        cell_sums = dataMatrix.sum(axis=1)
        to_keep = np.where(cell_sums > 0)[0]
        dataMatrix = dataMatrix.ix[dataMatrix.index[to_keep], :].astype(np.float32)

        #Remove empty genes
        print('Removing empty genes')
        gene_sums = dataMatrix.sum(axis=0)
        to_keep = np.where(gene_sums > 0)[0]
        dataMatrix = dataMatrix.ix[:, to_keep].astype(np.float32)

        # Construct class object
        scdata = cls( dataMatrix, data_type='sc-seq' )

        return scdata

    def filter_scseq_data(self, filter_cell_min=0, filter_cell_max=0, filter_gene_nonzero=None, filter_gene_mols=None):

        if filter_cell_min != filter_cell_max:
            sums = self.data.sum(axis=1)
            to_keep = np.intersect1d(np.where(sums >= filter_cell_min)[0], 
                                     np.where(sums <= filter_cell_max)[0])
            self.data = self.data.ix[self.data.index[to_keep], :].astype(np.float32)

        if filter_gene_nonzero != None:
            nonzero = self.data.astype(bool).sum(axis=0)
            to_keep = np.where(nonzero >= filter_gene_nonzero)[0]
            self.data = self.data.ix[:, to_keep].astype(np.float32)

        if filter_gene_mols != None:
            sums = self.data.sum(axis=0)
            to_keep = np.where(sums >= filter_gene_mols)[0]
            self.data = self.data.ix[:, to_keep].astype(np.float32)


    def normalize_scseq_data(self):
        """
        Normalize single cell RNA-seq data: Divide each cell by its molecule count 
        and multiply counts of cells by the median of the molecule counts
        :return: SCData
        """

        molecule_counts = self.data.sum(axis=1)
        data = self.data.div(molecule_counts, axis=0)\
            .mul(np.median(molecule_counts), axis=0)
        scdata = SCData(data=data, metadata=self.metadata)
        scdata._normalized = True

        # check that none of the genes are empty; if so remove them
        nonzero_genes = scdata.data.sum(axis=0) != 0
        scdata.data = scdata.data.ix[:, nonzero_genes].astype(np.float32)

        # set unnormalized_cell_sums
        self.library_sizes = molecule_counts
        scdata._library_sizes = molecule_counts

        return scdata


    def log_transform_scseq_data(self, pseudocount=0.1):
        self.data = np.log(np.add(self.data, pseudocount))
        
    def plot_molecules_per_cell_and_gene(self, fig=None, ax=None):

        height = 4
        width = 12
        fig = plt.figure(figsize=[width, height])
        gs = plt.GridSpec(1, 3)
        colsum = np.log10(self.data.sum(axis=0))
        rowsum = np.log10(self.data.sum(axis=1))
        for i in range(3):
            ax = plt.subplot(gs[0, i])

            if i == 0:
                print(np.min(rowsum))
                print(np.max(rowsum))
                n, bins, patches = ax.hist(rowsum, bins='auto')
                plt.xlabel('Molecules per cell (log10 scale)')
            elif i == 1:
                temp = np.log10(self.data.astype(bool).sum(axis=0))
                n, bins, patches = ax.hist(temp, bins='auto')
                plt.xlabel('Nonzero cells per gene (log10 scale)')
            else:
                n, bins, patches = ax.hist(colsum, bins='auto') 
                plt.xlabel('Molecules per gene (log10 scale)')
            plt.ylabel('Frequency')
            plt.tight_layout()
            ax.tick_params(axis='x', labelsize=8)

        return fig, ax


    def run_pca(self, n_components=100, random=True):
        """
        Principal component analysis of the data. 
        :param n_components: Number of components to project the data
        """
        
        solver = 'randomized'
        if random != True:
            solver = 'full'

        pca = PCA(n_components=n_components, svd_solver=solver)
        self.pca = pd.DataFrame(data=pca.fit_transform(self.data.values), index=self.data.index,
                                columns=['PC' + str(i) for i in range(1, n_components+1)])


    def plot_pca_variance_explained(self, n_components=30,
            fig=None, ax=None, ylim=(0, 100), random=True):
        """ Plot the variance explained by different principal components
        :param n_components: Number of components to show the variance
        :param ylim: y-axis limits
        :param fig: matplotlib Figure object
        :param ax: matplotlib Axis object
        :return: fig, ax
        """

        solver = 'randomized'
        if random != True:
            solver = 'full'
        pca = PCA(n_components=n_components, svd_solver=solver)
        pca.fit(self.data.values)

        fig, ax = get_fig(fig=fig, ax=ax)
        plt.plot(np.multiply(np.cumsum(pca.explained_variance_ratio_), 100))
        plt.ylim(ylim)
        plt.xlim((0, n_components))
        plt.xlabel('Components')
        plt.ylabel('Percent Variance explained')
        plt.title('Principal components')
        return fig, ax


    def run_tsne(self, n_components=50, perplexity=30, n_iter=1000, theta=0.5):
        """ Run tSNE on the data. tSNE is run on the principal component projections
        for single cell RNA-seq data and on the expression matrix for mass cytometry data
        :param n_components: Number of components to use for running tSNE for single cell 
        RNA-seq data. Ignored for mass cytometry
        :return: None
        """

        # Work on PCA projections if data is single cell RNA-seq
        if self.data_type == 'sc-seq':
            if self.pca is None:
                self.run_pca(n_components=n_components)
                data = self.pca
            elif n_components <= self.pca.shape[1]:
                data = self.pca.iloc[:, :n_components]
            else:
                pca = PCA(n_components=n_components, svd_solver='randomized')
                data = pca.fit_transform(self.data.values)
        else:
            data = self.data

        # Reduce perplexity if necessary
        perplexity_limit = 15
        if data.shape[0] < 100 and perplexity > perplexity_limit:
            print('Reducing perplexity to %d since there are <100 cells in the dataset. ' % perplexity_limit)
        tsne = TSNE(n_components=2, perplexity=perplexity, init='random', random_state=sum(data.shape), n_iter=n_iter, angle=theta) 
        self.tsne = pd.DataFrame(tsne.fit_transform(data),                       
								 index=self.data.index, columns=['tSNE1', 'tSNE2'])


    def plot_tsne(self, fig=None, ax=None, density=False, color=None, title='tSNE projection'):
        """Plot tSNE projections of the data
        :param fig: matplotlib Figure object
        :param ax: matplotlib Axis object
        :param title: Title for the plot
        """
        if self.tsne is None:
            raise RuntimeError('Please run tSNE using run_tsne before plotting ')
        fig, ax = get_fig(fig=fig, ax=ax)
        if isinstance(color, pd.Series):
            plt.scatter(self.tsne['tSNE1'], self.tsne['tSNE2'], s=size, 
                        c=color.values, edgecolors='none')
            plt.colorbar()
        elif density == True:
            # Calculate the point density
            xy = np.vstack([self.tsne['tSNE1'], self.tsne['tSNE2']])
            z = gaussian_kde(xy)(xy)

            # Sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            x, y, z = self.tsne['tSNE1'][idx], self.tsne['tSNE2'][idx], z[idx]

            plt.scatter(x, y, s=size, c=z, edgecolors='none')
            plt.colorbar()
        else:
            plt.scatter(self.tsne['tSNE1'], self.tsne['tSNE2'], s=size, edgecolors='none',
                        color=qualitative_colors(2)[1] if color == None else color)
        
        ax.set_title(title)
        plt.axis('tight')
        plt.tight_layout()
        return fig, ax


    def plot_tsne_by_cell_sizes(self, fig=None, ax=None, vmin=None, vmax=None):
        """Plot tSNE projections of the data with cells colored by molecule counts
        :param fig: matplotlib Figure object
        :param ax: matplotlib Axis object
        :param vmin: Minimum molecule count for plotting 
        :param vmax: Maximum molecule count for plotting 
        :param title: Title for the plot
        """
        if self.data_type == 'masscyt':
            raise RuntimeError( 'plot_tsne_by_cell_sizes is not applicable \n\
                for mass cytometry data. ' )

        fig, ax = get_fig(fig, ax)
        if self.tsne is None:
            raise RuntimeError('Please run run_tsne() before plotting.')
        if self._normalized:
            sizes = self.library_sizes
        else:
            sizes = self.data.sum(axis=1)
        plt.scatter(self.tsne['tSNE1'], self.tsne['tSNE2'], s=size, c=sizes, edgecolors='none')
        plt.colorbar()
        plt.axis('tight')
        plt.tight_layout()
        return fig, ax


    def run_diffusion_map(self, k=10, epsilon=1, distance_metric='euclidean',
        n_diffusion_components=10, n_pca_components=15, markers=None, ka=0, random_pca=True):
        """ Run diffusion maps on the data. Run on the principal component projections
        for single cell RNA-seq data and on the expression matrix for mass cytometry data
        :param k: Number of neighbors for graph construction to determine distances between cells
        :param epsilon: Gaussian standard deviation for converting distances to affinities
        :param n_diffusion_components: Number of diffusion components to Generalte
        :param n_pca_components: Number of components to use for running tSNE for single cell 
        RNA-seq data. Ignored for mass cytometry
        :return: None
        """

        if n_pca_components != None:
            if isinstance(self.pca, pd.DataFrame) and n_pca_components == len(self.pca.columns.values):
                data = self.pca
            else:
                data = magic.MAGIC.run_pca(self.data, n_components=n_pca_components, random=random_pca)
        else:
            data = self.data

        N = data.shape[0]

        # Nearest neighbors
        nbrs = NearestNeighbors(n_neighbors=k, metric=distance_metric).fit(data)
        distances, indices = nbrs.kneighbors(data)

        if ka > 0:
            print('Autotuning distances')
            for j in reversed(range(N)):
                temp = sorted(distances[j])
                lMaxTempIdxs = min(ka, len(temp))
                if lMaxTempIdxs == 0 or temp[lMaxTempIdxs] == 0:
                    distances[j] = 0
                else:
                    distances[j] = np.divide(distances[j], temp[lMaxTempIdxs])

        # Adjacency matrix
        rows = np.zeros(N * k, dtype=np.int32)
        cols = np.zeros(N * k, dtype=np.int32)
        dists = np.zeros(N * k)
        location = 0
        for i in range(N):
            inds = range(location, location + k)
            rows[inds] = indices[i, :]
            cols[inds] = i
            dists[inds] = distances[i, :]
            location += k
        if epsilon > 0:
            W = csr_matrix( (dists, (rows, cols)), shape=[N, N] )
        else:
            W = csr_matrix( (np.ones(dists.shape), (rows, cols)), shape=[N, N] )

        # Symmetrize W
        W = W + W.T

        if epsilon > 0:
            # Convert to affinity (with selfloops)
            rows, cols, dists = find(W)
            rows = np.append(rows, range(N))
            cols = np.append(cols, range(N))
            dists = np.append(dists/(epsilon ** 2), np.zeros(N))
            W = csr_matrix( (np.exp(-dists), (rows, cols)), shape=[N, N] )

        # Create D
        D = np.ravel(W.sum(axis = 1))
        D[D!=0] = 1/D[D!=0]

        #markov normalization
        T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(W)
        
        # Eigen value decomposition
        D, V = eigs(T, n_diffusion_components, tol=1e-4, maxiter=1000)
        D = np.real(D)
        V = np.real(V)
        inds = np.argsort(D)[::-1]
        D = D[inds]
        V = V[:, inds]
        V = T.dot(V)

        # Normalize
        for i in range(V.shape[1]):
            V[:, i] = V[:, i] / norm(V[:, i])
        V = np.round(V, 10)

        # Update object
        self.diffusion_eigenvectors = pd.DataFrame(V, index=self.data.index,
                                                   columns=['DC'+str(i) for i in range(n_diffusion_components)])
        self.diffusion_eigenvalues = pd.DataFrame(D)


    def plot_diffusion_components(self, other_data=None, title='Diffusion Components'):
        """ Plots the diffusion components on tSNE maps
        :return: fig, ax
        """
        if self.tsne is None:
            raise RuntimeError('Please run tSNE before plotting diffusion components.')
        if self.diffusion_eigenvectors is None:
            raise RuntimeError('Please run diffusion maps using run_diffusion_map before plotting')

        height = int(2 * np.ceil(self.diffusion_eigenvalues.shape[0] / 5))
        width = 10
        n_rows = int(height / 2)
        n_cols = int(width / 2)
        if other_data:
            height = height * 2
            n_rows = n_rows * 2
        fig = plt.figure(figsize=[width, height])
        gs = plt.GridSpec(n_rows, n_cols)

        for i in range(self.diffusion_eigenvectors.shape[1]):
            ax = plt.subplot(gs[i // n_cols, i % n_cols])

            plt.scatter(self.tsne['tSNE1'], self.tsne['tSNE2'], c=self.diffusion_eigenvectors[i],
                        edgecolors='none', s=size)

            plt.title( 'Component %d' % i, fontsize=10 )

        if other_data:
            for i in range(other_data.diffusion_eigenvectors.shape[1]):
                ax = plt.subplot(gs[(self.diffusion_eigenvectors.shape[1] + i) // n_cols, (self.diffusion_eigenvectors.shape[1] + i) % n_cols])

                plt.scatter(other_data.tsne['tSNE1'], other_data.tsne['tSNE2'], c=other_data.diffusion_eigenvectors[i],
                            edgecolors='none', s=size)

                plt.title( 'Component %d' % i, fontsize=10 )

        gs.tight_layout(fig)
        # fig.suptitle(title, fontsize=12)
        return fig, ax


    def plot_diffusion_eigen_vectors(self, fig=None, ax=None, title='Diffusion eigen vectors'):
        """ Plots the eigen values associated with diffusion components
        :return: fig, ax
        """
        if self.diffusion_eigenvectors is None:
            raise RuntimeError('Please run diffusion maps using run_diffusion_map before plotting')

        fig, ax = get_fig(fig=fig, ax=ax)
        ax.plot(np.ravel(self.diffusion_eigenvalues.values))
        plt.scatter( range(len(self.diffusion_eigenvalues)), 
            self._diffusion_eigenvalues, s=20, edgecolors='none', color='red' )
        plt.xlabel( 'Diffusion components')
        plt.ylabel('Eigen values')
        plt.title( title )
        plt.xlim([ -0.1, len(self.diffusion_eigenvalues) - 0.9])
        sns.despine(ax=ax)
        return fig, ax


    @staticmethod
    def _correlation(x: np.array, vals: np.array):
        x = x[:, np.newaxis]
        mu_x = x.mean()  # cells
        mu_vals = vals.mean(axis=0)  # cells by gene --> cells by genes
        sigma_x = x.std()
        sigma_vals = vals.std(axis=0)

        return ((vals * x).mean(axis=0) - mu_vals * mu_x) / (sigma_vals * sigma_x)


    def run_diffusion_map_correlations(self, components=None, no_cells=10):
        """ Determine gene expression correlations along diffusion components
        :param components: List of components to generate the correlations. All the components
        are used by default.
        :param no_cells: Window size for smoothing
        :return: None
        """
        if self.data_type == 'masscyt':
            raise RuntimeError('This function is designed to work for single cell RNA-seq')
        if self.diffusion_eigenvectors is None:
            raise RuntimeError('Please run diffusion maps using run_diffusion_map before determining correlations')

        if components is None:
            components = np.arange(self.diffusion_eigenvectors.shape[1])
        else:
            components = np.array(components)
        components = components[components != 0]

        # Container
        diffusion_map_correlations = np.empty((self.data.shape[1],
                                               self.diffusion_eigenvectors.shape[1]),
                                               dtype=np.float)
        for component_index in components:
            component_data = self.diffusion_eigenvectors.ix[:, component_index]

            order = self.data.index[np.argsort(component_data)]
            x = component_data[order].rolling(no_cells).mean()[no_cells:]
            # x = pd.rolling_mean(component_data[order], no_cells)[no_cells:]

            # this fancy indexing will copy self.data
            vals = self.data.ix[order, :].rolling(no_cells).mean()[no_cells:].values
            # vals = pd.rolling_mean(self.data.ix[order, :], no_cells, axis=0)[no_cells:]
            cor_res = self._correlation(x, vals)
            # assert cor_res.shape == (gene_shape,)
            diffusion_map_correlations[:, component_index] = self._correlation(x, vals)

        # this is sorted by order, need it in original order (reverse the sort)
        self.diffusion_map_correlations = pd.DataFrame(diffusion_map_correlations[:, components],
                            index=self.data.columns, columns=components)


    def plot_gene_component_correlations(
            self, components=None, fig=None, ax=None,
            title='Gene vs. Diffusion Component Correlations'):
        """ plots gene-component correlations for a subset of components

        :param components: Iterable of integer component numbers
        :param fig: Figure
        :param ax: Axis
        :param title: str, title for the plot
        :return: fig, ax
        """
        fig, ax = get_fig(fig=fig, ax=ax)
        if self.diffusion_map_correlations is None:
            raise RuntimeError('Please run determine_gene_diffusion_correlations() '
                               'before attempting to visualize the correlations.')

        if components is None:
            components = self.diffusion_map_correlations.columns
        colors = qualitative_colors(len(components))

        for c,color in zip(components, colors):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')  # catch experimental ipython widget warning
                sns.kdeplot(self.diffusion_map_correlations[c].fillna(0), label=c,
                            ax=ax, color=color)
        sns.despine(ax=ax)
        ax.set_title(title)
        ax.set_xlabel('correlation')
        ax.set_ylabel('gene density')
        plt.legend()
        return fig, ax


    # todo add option to plot phenograph cluster that these are being DE in.
    def plot_gene_expression(self, genes, other_data=None):
        """ Plot gene expression on tSNE maps
        :param genes: Iterable of strings to plot on tSNE        
        """
        if not isinstance(genes, dict):
            not_in_dataframe = set(genes).difference(self.extended_data.columns.get_level_values(1))
            if not_in_dataframe:
                if len(not_in_dataframe) < len(genes):
                    print('The following genes were either not observed in the experiment, '
                          'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
                else:
                    print('None of the listed genes were observed in the experiment, or the '
                          'wrong symbols were used.')
                return

        height = int(5 * np.ceil(len(genes) / 2))
        width = 10
        n_rows = int(height / 5)
        n_cols = 2
        if other_data:
            fig = plt.figure(figsize=[width, 2*(height+0.25)])
            gs = plt.GridSpec(2*n_rows, n_cols)
        else:
            fig = plt.figure(figsize=[width, height+0.25])
            gs = plt.GridSpec(n_rows, n_cols)

        axes = []
        for i, g in enumerate(genes):
            if not isinstance(genes, dict):
                g = self.extended_data.columns.values[np.where([g in col for col in self.extended_data.columns.values])[0]][0]
            ax = plt.subplot(gs[i // n_cols, i % n_cols])
            axes.append(ax)
            if isinstance(genes, dict):
                color = genes[g]
            else:
                color = self.extended_data[g]
            plt.scatter(self.tsne['tSNE1'], self.tsne['tSNE2'], c=color,
                        edgecolors='none', s=size)
            if isinstance(genes, dict):
                ax.set_title(g)
            else:           
                ax.set_title(g[1])
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            ax.set_xlabel('tSNE1')
            ax.set_ylabel('tSNE2')
            plt.colorbar()
            plt.axis('tight')

        if other_data:
            for i, g in enumerate(genes):
                if not isinstance(genes, dict):
                    g = self.extended_data.columns.values[np.where([g in col for col in self.extended_data.columns.values])[0]][0]
                
                ax = plt.subplot(gs[(n_rows*n_cols +i) // n_cols, (n_rows*n_cols +i) % n_cols])
                axes.append(ax)
                if isinstance(genes, dict):
                    color = genes[g]
                else:
                    color = other_data.extended_data[g]
                plt.scatter(other_data.tsne['tSNE1'], other_data.tsne['tSNE2'], c=color,
                            edgecolors='none', s=size)
  
                if isinstance(genes, dict):
                    ax.set_title(g)
                else:             
                    ax.set_title(g[1])
                ax.xaxis.set_major_locator(plt.NullLocator())
                ax.yaxis.set_major_locator(plt.NullLocator())
                ax.set_xlabel('tSNE1')
                ax.set_ylabel('tSNE2')
                plt.colorbar()
                plt.axis('tight')
        
        plt.tight_layout()
        return fig, axes


    def scatter_gene_expression(self, genes, density=False, color=None, fig=None, ax=None):
        """ 2D or 3D scatter plot of expression of selected genes
        :param genes: Iterable of strings to scatter
        """

        not_in_dataframe = set(genes).difference(self.extended_data.columns.get_level_values(1))
        if not_in_dataframe:
            if len(not_in_dataframe) < len(genes):
                print('The following genes were either not observed in the experiment, '
                      'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
            else:
                print('None of the listed genes were observed in the experiment, or the '
                      'wrong symbols were used.')
            return

        if len(genes) < 2 or len(genes) > 3:
            raise RuntimeError('Please specify either 2 or 3 genes to scatter.')

        for i in range(len(genes)):
            genes[i] = self.extended_data.columns.values[np.where([genes[i] in col for col in self.extended_data.columns.values])[0]][0]

        gui_3d_flag = True
        if ax == None:
            gui_3d_flag = False

        fig, ax = get_fig(fig=fig, ax=ax)
        if len(genes) == 2:
            if density == True:
                # Calculate the point density
                xy = np.vstack([self.extended_data[genes[0]], self.extended_data[genes[1]]])
                z = gaussian_kde(xy)(xy)

                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = self.extended_data[genes[0]][idx], self.extended_data[genes[1]][idx], z[idx]

                plt.scatter(x, y, s=size, c=z, edgecolors='none')
                ax.set_title('Color = density')
                plt.colorbar()
            elif isinstance(color, pd.Series):
                plt.scatter(self.extended_data[genes[0]], self.extended_data[genes[1]],
                            s=size, c=color, edgecolors='none')
                ax.set_title('Color = ' + color.name)
                plt.colorbar()
            elif color in self.extended_data.columns.get_level_values(1):
                color = self.extended_data.columns.values[np.where([color in col for col in self.extended_data.columns.values])[0]][0]
                plt.scatter(self.extended_data[genes[0]], self.extended_data[genes[1]],
                            s=size, c=self.extended_data[color], edgecolors='none')
                ax.set_title('Color = ' + color[1])
                plt.colorbar()
            else:
                plt.scatter(self.extended_data[genes[0]], self.extended_data[genes[1]], edgecolors='none',
                            s=size, color=qualitative_colors(2)[1] if color == None else color)
            ax.set_xlabel(genes[0][1])
            ax.set_ylabel(genes[1][1])

        else:
            if not gui_3d_flag:
                ax = fig.add_subplot(111, projection='3d')

            if density == True:
                xyz = np.vstack([self.extended_data[genes[0]],self.extended_data[genes[1]],
                                 self.extended_data[genes[2]]])
                kde = gaussian_kde(xyz)
                density = kde(xyz)

                p = ax.scatter(self.extended_data[genes[0]], self.extended_data[genes[1]], self.extended_data[genes[2]],
                           s=size, c=density, edgecolors='none')
                ax.set_title('Color = density')
                fig.colorbar(p)
            elif isinstance(color, pd.Series):
                p = ax.scatter(self.extended_data[genes[0]], self.extended_data[genes[1]],
                           self.extended_data[genes[2]], s=size, c=color, edgecolors='none')
                ax.set_title('Color = ' + color.name)
                fig.colorbar(p)
            elif color in self.extended_data.columns.get_level_values(1):
                color = self.extended_data.columns.values[np.where([color in col for col in self.extended_data.columns.values])[0]][0]
                p = ax.scatter(self.extended_data[genes[0]], self.extended_data[genes[1]],
                           self.extended_data[genes[2]], s=size, c=self.extended_data[color], edgecolors='none')
                ax.set_title('Color = ' + color[1])
                fig.colorbar(p)
            else:
                p = ax.scatter(self.extended_data[genes[0]], self.extended_data[genes[1]], self.extended_data[genes[2]], 
                           edgecolors='none', s=size, color=qualitative_colors(2)[1] if color == None else color)
            ax.set_xlabel(genes[0][1])
            ax.set_ylabel(genes[1][1])
            ax.set_zlabel(genes[2][1])
            ax.view_init(15,55)
        
        plt.axis('tight')
        plt.tight_layout()
        return fig, ax


    def scatter_gene_expression_against_other_data(self, genes, other_data, density=False, color=None, fig=None, ax=None):

        not_in_dataframe = set(genes).difference(self.extended_data.columns.get_level_values(1))
        if not_in_dataframe:
            if len(not_in_dataframe) < len(genes):
                print('The following genes were either not observed in the experiment, '
                      'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
            else:
                print('None of the listed genes were observed in the experiment, or the '
                      'wrong symbols were used.')
            return

        height = int(4 * np.ceil(len(genes) / 3))
        width = 12 if len(genes) >= 3 else 4*len(genes)
        n_rows = int(height / 4)
        n_cols = int(width / 4)
        fig = plt.figure(figsize=[width, height])
        gs = plt.GridSpec(n_rows, n_cols)

        axes = []
        for i, g in enumerate(genes):
            g = self.extended_data.columns.values[np.where([g in col for col in self.extended_data.columns.values])[0]][0]
            ax = plt.subplot(gs[i // n_cols, i % n_cols])
            axes.append(ax)
            if density == True:
                # Calculate the point density
                xy = np.vstack([self.extended_data[g], other_data.extended_data[g]])
                z = gaussian_kde(xy)(xy)

                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = self.extended_data[g][idx], other_data.extended_data[g][idx], z[idx]

                plt.scatter(x, y, s=size, c=z, edgecolors='none')
                plt.colorbar()
            elif isinstance(color, pd.Series):
                plt.scatter(self.extended_data[g], other_data.extended_data[g], s=size, c=color, edgecolors='none') 
                plt.colorbar()
            elif color in self.extended_data.columns.get_level_values(1):
                color = self.extended_data.columns.values[np.where([color in col for col in self.extended_data.columns.values])[0]][0]
                plt.scatter(self.extended_data[g], other_data.extended_data[g], s=size, c=self.extended_data[color], edgecolors='none')
                plt.colorbar()
            else:
                plt.scatter(self.extended_data[g], other_data.extended_data[g], s=size, edgecolors='none',
                            color=qualitative_colors(2)[1] if color == None else color)  
            plt.axis('tight')
                       
        gs.tight_layout(fig, pad=3, h_pad=3, w_pad=3)

        return fig, axes


    def run_magic(self, n_pca_components=20, random_pca=True, t=6, k=30, ka=10, epsilon=1, rescale_percent=99):

        new_data = magic.MAGIC.magic(self.data.values, n_pca_components=n_pca_components, random_pca=random_pca, t=t, 
                                     k=k, ka=ka, epsilon=epsilon, rescale=rescale_percent)

        new_data = pd.DataFrame(new_data, index=self.data.index, columns=self.data.columns)

        # Construct class object
        scdata = magic.mg.SCData(new_data, data_type=self.data_type)
        self.magic = scdata

    def concatenate_data(self, other_data_sets, join='outer', axis=0, names=[]):

        #concatenate dataframes
        temp = self.data.copy()
        if axis == 0:
            temp.index = [names[0] + ' ' + i for i in self.data.index]
        else:
            temp.columns = [names[0] + ' ' + i for i in self.data.columns]
        dfs = [temp]
        count = 0
        for data_set in other_data_sets:
            count += 1
            temp = data_set.data.copy()
            if axis == 0:
                temp.index = [names[count] + ' ' + i for i in data_set.data.index]
            else:
                temp.columns = [names[count] + ' ' + i for i in self.data.columns]
            dfs.append(temp)
        df_concat = pd.concat(dfs, join=join, axis=axis)

        scdata = magic.mg.SCData(df_concat)
        return scdata
