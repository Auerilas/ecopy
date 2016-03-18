import numpy as np
from pandas import DataFrame

def div_partition(x, method='shannon', breakNA=True, weights=None):
	'''
	Docstring for function ecopy.diversity
	========================
	Decomposes diversity measures into alpha, beta, and gamma components.

	Use
	----
	div_partition(x, method='shannon', breakNA=True)

	Returns a tuple of alpha, beta, and gamma diversity

	Parameters
	----------
	x:  numpy array or pandas dataframe with observations as rows
		and descriptors as columns
	method: a method used for calculating species diversity. See the 
		diversity() function for details. Can only be 'shannon', 
		'gini-simpson', 'simpson', or 'spRich'.
	breakNA: should the process halt if the matrix contains any NAs?
		if False, then NA's undergo pairwise deletion during distance calculation,
		such that when calculating the distance between two rows, if any
		species is missing from a row, that species is removed from both rows
	weights: weights given for each row (site). Defaults to the 
		sum of each row divided by the sum of the matrix. This yields
		weights based on the number of individuals in a site for raw
		abundance data or equal weights for relative abundance data.

	Example
	--------
	import ecopy as ep
	varespec = ep.load_data('varespec')
	
	D_alpha, D_beta, D_gamma = ep.div_partition(varespec, 'shannon')
	'''
	listofmethods = ['shannon', 'gini-simpson', 'simpson', 'dominance', 'spRich', 'even']
	if not isinstance(breakNA, bool):
		msg = 'removaNA argument must be boolean'
		raise ValueError(msg)
	if method not in listofmethods:
		msg = 'method argument {0!s} is not an accepted metric'.format(method)
		raise ValueError(msg)
	if not isinstance(x, (DataFrame, np.ndarray)):
		msg = 'x argument must be a numpy array or pandas dataframe'
		raise ValueError(msg)
	if isinstance(x, DataFrame):
		if (x.dtypes == 'object').any():
			msg = 'DataFrame can only contain numeric values'
		if breakNA:
			if x.isnull().any().any():
				msg = 'DataFrame contains null values'
				raise ValueError(msg)
		if (x<0).any().any():
			msg = 'DataFrame contains negative values'
			raise ValueError(msg)
		z = np.array(x, 'float')
	if isinstance(x, np.ndarray):
		if breakNA:
			if np.isnan(np.sum(x)):
				msg = 'Array contains null values'
				raise ValueError(msg)
		if (x < 0).any():
			msg = 'Array contains negative values'
			raise ValueError(msg)
		z = np.array(x, 'float')
	if weights is None:
		weights = z.sum(axis=1)/z.sum()
	w = np.array(weights)/np.array(weights).sum()
	relMat = z / z.sum(axis=1)[:,np.newaxis]
	totalList = z.sum(axis=0)
	totalRel = totalList / totalList.sum()
	if method=='shannon':
		div = np.apply_along_axis(shannonFunc, 1, relMat)
		H_alpha = (w*div).sum() / w.sum()
		H_gamma = shannonFunc(totalRel)
		D_alpha = np.exp(H_alpha)
		D_gamma = np.exp(H_gamma)
		D_beta = D_gamma / D_alpha
		return D_alpha, D_beta, D_gamma
	if method=='gini-simpson':
		div = np.apply_along_axis(giniFunc, 1, relMat)
		H_alpha = (w*div).sum()/w.sum()
		H_gamma = giniFunc(totalRel)
		D_alpha = 1./(1.-H_alpha)
		D_gamma = 1./(1.-H_gamma)
		D_beta = D_gamma / D_alpha
		return D_alpha, D_beta, D_gamma
	if method=='simpson':
		div = np.apply_along_axis(simpson, 1, relMat)
		H_alpha = (w*div).sum()/w.sum()
		H_gamma = simpson(totalRel)
		D_alpha = 1./H_alpha
		D_gamma = 1./H_gamma
		D_beta = D_gamma / D_alpha
		return D_alpha, D_beta, D_gamma
	if method=='spRich':
		div = np.apply_along_axis(richness, 1, relMat)
		H_alpha = (w*div).sum()/w.sum()
		H_gamma = richness(totalRel)
		D_alpha = H_alpha
		D_gamma = H_gamma
		D_beta = D_gamma / D_alpha
		return D_alpha, D_beta, D_gamma

def shannonFunc(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	t = t[t!=0]
	H = -np.sum( t*np.log(t) )
	return H

def giniFunc(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	D = 1 - np.sum( t**2 )
	return D

def simpson(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	D = np.sum( t**2 )
	return D

def richness(y):
	notabs = ~np.isnan(y)
	t = y[notabs]
	D = np.sum(t!=0)
	return float(D)
