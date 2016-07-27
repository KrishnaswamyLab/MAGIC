import numpy as np
import pandas as pd
from scipy import sparse
import wishbone

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
	data_new = np.array(np.dot(data, L_t))

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
    
	return data_new, L_t

def run_magic(data, n_pca_components=None, t=8, knn=20, epsilon=1, rescale=True):

	if n_pca_components != None:
		#project data onto n pca components
		if isinstance(data, wishbone.wb.SCData):
			data = data.data
		if isinstance(data, pd.DataFrame):
			pca_projected_data = run_pca(data.values, n_pca_components)
		else:
			pca_projected_data = run_pca(data, n_pca_components)
	else:
		pca_projected_data = data

	#run diffusion maps to get markov matrix
	diffusion_map = run_diffusion_map(pca_projected_data, knn=knn, normalization='markov', epsilon=epsilon, distance_metric='euclidean')

	#get imputed data matrix
	new_data, L_t = impute_fast(data, diffusion_map['T'], t, rescale_to_max=rescale)

	if isinstance(data, pd.DataFrame):
		new_data = pd.DataFrame(new_data, index=data.index, columns=data.columns)
	else:
		new_data = pd.DataFrame(new_data)

	return new_data, L_t


