from pandas import DataFrame
import numpy as np

class Mantel(object):
	'''
	Docstring for function ecopy.Mantel
	====================
	Conducts a Mantel test for association between two square, symmetric, non-negative
		distance matrices. Allows for partial Mantel tests if given a third conditioning 
		matrix

	Use
	----
	Mantel(d1, d2, d_condition=None, test='pearson', tail='both', nperm=999)

	Returns an object of class Mantel

	Parameters
	----------
	d1:  First distance matrix
	d2: Second distance matrix for comparison
	d_condition: Third distance matrix to condition correlations on (i.e. partial Mantel test)
	test: 'pearson' operates on standardized variables, 
		'spearman' operates on standardized ranks
 	tail: 'greater' tests the one-tailed hypothesis that correlation is 
 		greater than predicted. 'lower' tests hypothsis that correlation
 		is lower than predicted. 'both' is a two-tailed test
 	nperm: Number of permutations for the test.

	Attributes (see online documentation for descriptions)
	---------
	r_obs: Observed correlation statistic
	pval: p-value for the given hypothesis
	tail: The tested hypothesis
	test: Which of the statistics used, 'pearson' or 'spearman'
	perm: Number of permutations

	Methods
	--------
	summary(): provides a summary of test results

	Example
	--------
	import ecopy as ep

	v1 = ep.load_data('varespec')
	v2 = ep.load_data('varechem')
	v2 = v2.apply(lambda x: (x - x.mean())/x.std(), axis=0)

	dist1 = ep.distance(v1, 'bray')
	dist2 = ep.distance(v2, 'euclidean')

	mant = ep.Mantel(dist1, dist2)
	print(mant.summary())
	'''
	def __init__(self, d1, d2, d_condition = None, test='pearson', tail='both', nperm=999):
		if not isinstance(d1, (np.ndarray, DataFrame)):
			msg = 'Matrix d1 must be a numpy.ndarray or pandas.DataFrame'
		if not isinstance(d2, (np.ndarray, DataFrame)):
			msg = 'Matrix d2 must be a numpy.ndarray or pandas.DataFrame'
		if isinstance(d1, DataFrame):
			d1 = np.array(d1)
		if isinstance(d2, DataFrame):
			d2 = np.array(d1)
		if np.any(d2 < 0):
			msg ='Matrix d2 cannot have negative values'
			raise ValueError(msg)
		if d1.shape[0] != d1.shape[1]:
			msg = 'Matrix d1 must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if d2.shape[0] != d2.shape[1]:
			msg = 'Matrix d2 must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if not np.allclose(d1.T, d1):
			msg = 'Matrix d1 must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if not np.allclose(d2.T, d2):
			msg = 'Matrix d2 must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if d1.shape[0] != d2.shape[0]:
			msg = 'Matrices must have same dimensions'
			raise ValueError(msg)
		if test not in ['pearson', 'spearman']:
			msg = 'test must be pearson or spearman'
			raise ValueError(msg)
		if tail not in ['greater', 'lower', 'both']:
			msg = 'tail must be greater, lower, both'
			raise ValueError(msg)
		if nperm < 2:
			msg = 'nperm must be > 2'
			raise ValueError(msg)
		r_perm = np.empty(nperm)
		if d_condition is not None:
			if d_condition.shape[0] != d_condition.shape[1]:
				msg = 'Conditioning matrix must be a square, symmetric distance matrix'
				raise ValueError(msg)
			if not np.allclose(d_condition.T, d_condition):
				msg = 'Conditioning matrix must be a square, symmetric distance matrix'
				raise ValueError(msg)
			resMat = residCalc(d1, d_condition)
			partR = partialMantel(resMat, d2, d_condition)
			self.robs = partR[0,1] 
			for i in range(nperm):
				residMat_perm = permuteFunc(resMat)
				r_star = partialMantel(residMat_perm, d2, d_condition)
				r_perm[i] = r_star[0,1]
		if test is 'pearson':
			self.r_obs = manFunc_pears(d1, d2)
			for i in range(nperm):
				d1_perm = permuteFunc(d1)
				r_perm[i] = manFunc_pears(d1_perm, d2)
		if test is 'spearman':
			self.r_obs = manFunc_spear(d1, d2)
			for i in range(nperm):
				d1_perm = permuteFunc(d1)
				r_perm[i] = manFunc_spear(d1_perm, d2)
		if tail is 'greater':
			self.pval = np.mean(r_perm > self.r_obs)
		elif tail is 'lower':
			self.pval = np.mean(r_perm < self.r_obs)
		else:
			self.pval = np.mean(np.abs(r_perm) > self.r_obs)
		self.test = test
		self.tail = tail
		self.perm = nperm

	def summary(self):
		summ = '\n{0} Mantel Test\nHypothesis = {1}\n\nObserved r = {2:.3}\tp = {3:.3}\n{4} permutations'.format(self.test.title(), self.tail, self.r_obs, self.pval, self.perm)
		return summ



def manFunc_pears(x,y):
	n = x.shape[0]
	denom = n*(n-1)/2. - 1.
	ui = np.triu_indices(x.shape[0])
	x2 = x.copy()
	y2 = y.copy()
	x2[ui] = np.nan
	y2[ui] = np.nan
	x_flat = x2.ravel()
	y_flat = y2.ravel()
	x_flat = x_flat[~np.isnan(x_flat)]
	y_flat = y_flat[~np.isnan(y_flat)]
	x_flat = (x_flat - x_flat.mean())/x_flat.std(ddof=1)
	y_flat = (y_flat - y_flat.mean())/y_flat.std(ddof=1)
	r = x_flat.dot(y_flat) / denom
	return r

def manFunc_spear(x,y):
	n = x.shape[0]
	denom = n*(n-1)/2. - 1.
	ui = np.triu_indices(x.shape[0])
	x2 = x.copy()
	y2 = y.copy()
	x2[ui] = np.nan
	y2[ui] = np.nan
	x_flat = x2.ravel()
	y_flat = y2.ravel()
	x_flat = x_flat[~np.isnan(x_flat)]
	y_flat = y_flat[~np.isnan(y_flat)]
	x_flat = x_flat.argsort().argsort()
	y_flat = y_flat.argsort().argsort()
	x_flat = (x_flat - x_flat.mean())/x_flat.std(ddof=1)
	y_flat = (y_flat - y_flat.mean())/y_flat.std(ddof=1)
	r = x_flat.dot(y_flat) / denom
	return r

def permuteFunc(x):
	idx = np.random.choice(x.shape[0], x.shape[0])
	xperm = x[idx,:]
	for j in range(xperm.shape[0]):
		xperm[j, idx[j]] = xperm[j, j]
	np.fill_diagonal(xperm, 0)
	return xperm

def residCalc(y, z):
	y2 = y.copy()
	z2 = z.copy()
	ui = np.triu_indices(y2.shape[0])
	iy, ix = np.indices(y2.shape)
	y2[ui] = np.nan
	z2[ui] = np.nan
	y_flat = y2.ravel()
	z_flat = z2.ravel()
	ix_flat = ix.ravel()
	iy_flat = iy.ravel()
	ix_flat = ix_flat[~np.isnan(y_flat)]
	iy_flat = iy_flat[~np.isnan(y_flat)]
	y_flat = y_flat[~np.isnan(y_flat)]
	z_flat = z_flat[~np.isnan(z_flat)]
	Z = np.array([np.ones(len(z_flat)), z_flat]).T
	params = np.linalg.lstsq(Z, y_flat)[0]
	Res = y_flat - Z.dot(params)
	ResMat = np.zeros(y2.shape)
	for i in range(len(ix_flat)):
		ResMat[ix_flat[i], iy_flat[i]] = Res[i]
	for i in range(y2.shape[0]):
		for j in range(y2.shape[0]):
			ResMat[j,i] = ResMat[i,j]
	return ResMat

def partialMantel(x, y, z):
	R = np.eye(3)
	d = [x, y, z]
	for i in range(3):
		for j in range(3):
			R[i,j] = manFunc_pears(d[i], d[j])
	R11 = R[:2, :2]
	R12 = R[:2, 2].reshape(2, 1)
	R21 = R12.T
	Rpart = R11 - R12.dot(R21)
	Dpart = np.diag(np.diag(Rpart)**-0.5)
	R12_3 = Dpart.dot(Rpart).dot(Dpart)
	return R12_3
