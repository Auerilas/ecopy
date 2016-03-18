import numpy as np
from pandas import DataFrame
from scipy.optimize import minimize

def spatial_median(X):
	"""
	Docstring for function ecopy.spatial_median
	====================
	Calculates the spatial median of a multivariate dataset.

	Use
	----
	spatial_median(X)

	Returns a numpy.ndarray

	Parameters
	----------
	X: A numpy.ndarray or pandas.DataFrame

	Example
	--------	
	import ecopy as ep
	from scipy.stats import multivariate_normal

	np.random.seed(654321)
	cov = np.diag([3.,5.,2.])
	data = multivariate_normal.rvs([0,0,0], cov, (100,))
	spatialMed = ep.spatial_median(data)
	"""
	if not isinstance(X, (DataFrame, np.ndarray)):
		msg = 'X must be a pandas.DataFrame or numpy.ndarray'
		raise ValueError(msg)
	X = np.array(X)
	medians = np.apply_along_axis(np.median, 0, X)
	output = minimize(medianSearch, medians, args=(X,), method='BFGS')
	finalMedians = output['x']
	return(finalMedians)

def medianSearch(params, Z):
	med_iter = params
	eucDist = np.apply_along_axis(lambda z: np.sqrt(np.sum((z-med_iter)**2)), 1, Z)
	err = eucDist
	return(np.mean(err))