import numpy as np
import time
from scipy.sparse import csr_matrix

def bimarkov(K, max_iters=100, abs_error=0.00001):
	"""BIMARKOV computes the bimarkov normalization function p(x) for the
	nonnegative, symemtric  kernel K(x,y) using an iterative scheme.  The
	function p(x) is the unique function s.t.

    diag(1./sqrt(p)*K*diag(1./sqrt(p))

	is both row and column stochastic.  Note that a bimarkov kernel
	necessarily has L^2 norm 1."""

	start = time.process_time()

	if K.size == 0:
		return
		
	#process input
	if K.shape[0] != K.shape[1]:
		print('Bimarkov.py: kernel must be NxN\n')
		return

	N = K.shape[0]

	print('Bimarkov: normalizing %dx%d kernel ... ', N, N)

	# initialize
	p = np.ones(N)

	# iterative
	for i in range(max_iters):

		S = np.ravel(W.sum(axis = 1))
		err = np.max(np.absolute(1.0-np.max(S)), np.absolute(1.0-np.min(S)))

		print('iter %03d error: 1e%1.4f', i, np.log10(err))

		if err < abs_error:
			if np.log10(err) < 0:
				print('\b')
			print('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
			break

		D = csr_matrix((np.divide(1, np.sqrt(S)), (range(N), range(N))), shape=[N, N])
		p *= S 
		K = D.dot(K).dot(D)

		if np.log10(err) < 0:
			print('\b')
		print('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')

	print('\b done (%4.2f ms, iters = %d, err = %g) \n', (time.process_time()-start)*1000.0, i, err)

	# iron out numerical errors
	T = (K + K.T)/2

	return T, p



