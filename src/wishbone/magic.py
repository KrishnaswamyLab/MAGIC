import numpy as np
from scipy import sparse
from GraphDiffusion.graph_diffusion import run_diffusion_map

def run_pca(X, n_components):

	# Make sure data is zero mean
    X = np.subtract(X, np.amin(X))
    X = np.divide(X, np.amax(X))

    # Compute covariance matrix
    if (X.shape[1] < X.shape[0]):
        C = np.cov(X, rowvar=0)
    # if N>D, we better use this matrix for the eigendecomposition
    else:
        C = np.multiply((1/X.shape[0]), np.dot(X, X.T))

    # Perform eigendecomposition of C
    C[np.where(np.isnan(C))] = 0
    C[np.where(np.isinf(C))] = 0
    l, M = np.linalg.eig(C)

    # Sort eigenvectors in descending order
    ind = np.argsort(l)[::-1]
    l = l[ind]
    if n_components < 1:
        n_components = np.where(np.cumsum(np.divide(l, np.sum(l)), axis=0) >= n_components)[0][0] + 1
        print('Embedding into ' + str(n_components) + ' dimensions.')
    if n_components > M.shape[1]:
        n_components = M.shape[1]
        print('Target dimensionality reduced to ' + str(n_components) + '.')

    M = M[:, ind[:n_components]]
    l = l[:n_components]

    # Apply mapping on the data
    if X.shape[1] >= X.shape[0]:
        M = np.multiply(np.dot(X.T, M), (1 / np.sqrt(X.shape[0] * l)).T)

    loadings = pd.DataFrame(data=M, index=self.molecules.columns)
    l = pd.DataFrame(l)

    return loadings

def impute_fast(data, L, t, rescale_to_max, L_t=None, tprev=None):

	#convert L to full matrix
	L = L.todense()

	#L^t
	print('MAGIC: L_t = L^t')
	if L_t == None:
		L_t = np.linalg.matrix_power(L, t)
	else:
		L_t = np.dot(L_t, np.linalg.matrix_power(L, t-tprev))

	print('MAGIC: data_new = L_t * data')
	data_new = np.dot(L_t, data)

	#rescale data to 99th percentile
	if rescale_to_max == True:
		M99 = np.percentile(data, 99, axis=0)
		M100 = data.max(axis=0)
		indices = np.where(M99 == 0)[0]
		M99[indices] = M100[indices]
		M99_new = np.percentile(data_new, 99, axis=0)
		M100_new = data_new.max(axis=0)
		indices = np.where(M99_new == 0)[0]
		M99_new[indices] = M100_new[indices]
		max_ratio = np.divide(M99, M99_new)
		data_new = np.multiply(data_new, np.matlib.repmat(max_ratio, len(data), 1))

def run_magic(data, n_pca_components=20, t=8, knn=30, epsilon=0, rescale=True):

	#project data onto n pca components
	pca_projected_data = run_pca(data, n_pca_components)

	#run diffusion maps to get markov matrix
	diffusion_map = run_diffusion_map(pca_projected_data, knn=knn, normalization='markov', epsilon=epsilon, distance_metric='euclidean')

	#get imputed data matrix
	return impute_fast(data, diffusion_map['T'], t, rescale_to_max=rescale)
