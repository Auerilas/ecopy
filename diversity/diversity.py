import numpy as np
from pandas import DataFrame 

def diversity(x, method='shannon', breakNA=True):
	'''
	Docstring for function ecopy.diversity
	========================
	Computes a given diversity index for a site x species matrix.
		All indices computed along rows (axis = 1)

	Use
	----
	ecopy.diversity(x, method='shannon', breakNA=True)

	Returns a numpy.ndarray or a pandas.Series

	Parameters
	----------
	x:  numpy array or pandas dataframe with observations as rows
		and descriptors as columns
	method: a method used for calculating species diversity
		shannon: Shannon's H
			H = -sum( p_i * log(p_i) )
			where p_i is relative abundance of species i
		simpson: Simpson's D
			D = 1 - sum( p_i^2 )
		invSimpson: Inverse Simpson's D
			D = 1/sum(p_i^2)
		dominance: Dominance index
			D = max(p_i)
		spRich: Species richness
			D = Number of non-zero observations
		even: Evenness of a sample. This is Shannon's H divided by the
			log of species richness

	breakNA: should the process halt if the matrix contains any NAs?
		if False, then NA's undergo pairwise deletion during distance calculation,
		such that when calculating the distance between two rows, if any
		species is missing from a row, that species is removed from both rows

	Example
	--------
	import ecopy as ep
	varespec = ep.load_data('varespec'')
	
	div = ep.diversity(varespec, 'shannon')
	'''
	listofmethods = ['shannon', 'simpson', 'invSimpson', 'dominance', 'spRich', 'even']
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
		if method=='shannon':
			div = x.apply(shannonFunc, axis=1)
			return div
		if method=='simpson':
			div = x.apply(simpFunc, axis=1)
			return div
		if method=='invSimpson':
			div = x.apply(simpInv, axis=1)
			return div
		if method=='dominance':
			div = x.apply(dom, axis=1)
			return div
		if method=='spRich':
			div = x.apply(richness, axis=1)
			return div
		if method=='even':
			div = x.apply(evenFunc, axis=1)
			return div
	if isinstance(x, np.ndarray):
		if breakNA:
			if np.isnan(np.sum(x)):
				msg = 'Array contains null values'
				raise ValueError(msg)
		if (x < 0).any():
			msg = 'Array contains negative values'
			raise ValueError(msg)
		if method=='shannon':
			div = np.apply_along_axis(shannonFunc, 1, x)
			return div
		if method=='simpson':
			div = np.apply_along_axis(simpFunc, 1, x)
			return div
		if method=='invSimpson':
			div = np.apply_along_axis(simpInv, 1, x)
			return div
		if method=='dominance':
			div = np.apply_along_axis(dom, 1, x)
			return div
		if method=='spRich':
			div = np.apply_along_axis(richness, 1, x)
			return div
		if method=='even':
			div = np.apply_along_axis(evenFunc, 1, x)
			return div

def shannonFunc(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	t = t[t!=0]
	H = -np.sum( t*np.log(t) )
	return H

def simpFunc(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	D = 1 - np.sum( t**2 )
	return D

def simpInv(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	D = 1. / np.sum( t**2 )
	return D

def dom(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	D = np.max(t)
	return D

def richness(y):
	notabs = ~np.isnan(y)
	t = y[notabs]
	D = np.sum(t!=0)
	return float(D)

def evenFunc(y):
	notabs = ~np.isnan(y)
	t = y[notabs] / np.sum(y[notabs])
	n = float(np.sum(t!=0))
	t = t[t!=0]
	H = -np.sum( t*np.log(t) )
	return H/np.log(n)
