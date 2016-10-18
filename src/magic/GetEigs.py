import numpy as np
from scipy.sparse.linalg import eigs
from numpy.linalg import norm

def GetEigs(T, k, P, take_diagonal=0):
	""" return k largest magnitude eigenvalues for the matrix T.
	:param T: Matrix to find eigen values/vectors of
	:param k: number of eigen values/vectors to return
	:param P: in the case of symmetric normalizations, 
	this is the NxN diagonal matrix which relates the nonsymmetric 
	version to the symmetric form via conjugation
	:param take_diagonal: if 1, returns the eigenvalues as a vector rather than as a diagonal matrix.
	"""
	D, V = eigs(T, k, tol=1e-4, maxiter=1000)
	D = np.real(D)
	V = np.real(V)
	inds = np.argsort(D)[::-1]
	D = D[inds]
	V = V[:, inds]
	if P is not None:
	    V = P.dot(V)

	# Normalize
	for i in range(V.shape[1]):
	    V[:, i] = V[:, i] / norm(V[:, i])
	V = np.round(V, 10)

	if take_diagonal == 0:
		D = np.diag(D)
		
	return V, D
