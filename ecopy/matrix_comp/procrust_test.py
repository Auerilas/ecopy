import numpy as np
from pandas import DataFrame

class procrustes_test(object):
	"""
	Docstring for function ecopy.procrustes_test
	====================
	Conducts permutation procrustes test of relationship
		between two non-diagonal (raw) matrices

	Use
	----
	procrustes_test(mat1, mat2, nperm)

	Returns an object of class procrustes_test

	Parameters
	----------
	mat1: A raw site x species matrix (or any object x descriptor)
	mat2: A raw site x descriptor matrix (or any object x descriptor)
	nperm: Number of permutations

	Attributes (see online documentation for descriptions)
	---------
	m12_obs: Observed test statistic, m12**2
	pval: p-value
	perm: Number of permutations


	Methods
	--------
	summary(): provides a summary of test results

	Example
	--------
	import ecopy as ep

	d1 = ep.load_data('varespec')
	d2 = ep.load_data('varechem')

	d = ep.procrustes_test(d1, d2)
	print(d.summary())
	"""
	def __init__(self, mat1, mat2, nperm=999):
		if isinstance(mat1, DataFrame):
			X = np.array(mat1).astype('float')
		else:
			X = mat1.astype('float')
		if isinstance(mat2, DataFrame):
			Y = np.array(mat2).astype('float')
		else:
			Y = mat2.astype('float')
		if X.shape[0] != Y.shape[0]:
			msg = 'Matrices must have the same number of rows'
			raise ValueError(msg)
		X_cent = np.apply_along_axis(lambda x: x - x.mean(), 0, X)
		Y_cent = np.apply_along_axis(lambda y: y - y.mean(), 0, Y)
		X_cent = X_cent / np.sqrt(np.sum(X_cent**2))
		Y_cent = Y_cent / np.sqrt(np.sum(Y_cent**2))
		W = np.sum(np.linalg.svd(X_cent.T.dot(Y_cent), compute_uv=0))
		self.m12_obs = 1 - W**2
		m12_perm = np.zeros(nperm)
		i = 0
		while i < nperm:
			idx = np.random.permutation(range(X_cent.shape[0]))
			X_perm = X_cent[idx,:]
			W_perm = np.sum(np.linalg.svd(X_perm.T.dot(Y_cent), compute_uv=0))
			m12_perm[i] = 1 - W_perm**2
			i += 1
		self.pval = np.mean(m12_perm < self.m12_obs)
		self.perm = nperm

	def summary(self):
		summ = '\nm12 squared = {0:.3}\np = {1:.3}\npermutations = {2}'.format(self.m12_obs, self.pval, self.perm)
		return summ


