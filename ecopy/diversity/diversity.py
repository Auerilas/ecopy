import numpy as np
from pandas import DataFrame, Series

def diversity(x, method='shannon', breakNA=True, num_equiv=True):
	'''
	Docstring for function ecopy.diversity
	========================
	Computes a given diversity index for a site x species matrix.
		All indices computed along rows (axis = 1)

	Use
	----
	diversity(x, method='shannon', breakNA=True, num_equiv=True)

	Returns a numpy.ndarray or a pandas.Series

	Parameters
	----------
	x:  numpy array or pandas dataframe with observations as rows
		and descriptors as columns
	method: a method used for calculating species diversity
		shannon: Shannon's H
			H = -sum( p_i * log(p_i) )
			where p_i is relative abundance of species i
		gini-simpson: Gini-Simpson
			D = 1 - sum( p_i^2 )
		simpson: Simpson's D
			D = sum(p_i^2)
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
	num_equiv: Whether or not the diversity index return species number equivalents,
		i.e. the number of species of identical abundance. This has better properties 
		than raw diversity. The number equivalents are as follows:

		shannon: np.exp(H)
		gini-simpson: 1 / (1-D)
		simpson: 1/D

	Example
	--------
	import ecopy as ep
	varespec = ep.load_data('varespec')
	
	div = ep.diversity(varespec, 'shannon')
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
	z = z / z.sum(axis=1)[:,np.newaxis]
	if method=='shannon':
		div = np.apply_along_axis(shannonFunc, 1, z)
		if num_equiv:
			div = np.exp(div)
		return div
	if method=='gini-simpson':
		div = np.apply_along_axis(giniFunc, 1, z)
		if num_equiv:
			div = 1./(1.-div)
		return div
	if method=='simpson':
		div = np.apply_along_axis(simpson, 1, z)
		if num_equiv:
			div = 1./div
		return div
	if method=='dominance':
		div = np.apply_along_axis(dom, 1, z)
		return div
	if method=='spRich':
		div = np.apply_along_axis(richness, 1, z)
		return div
	if method=='even':
		div = np.apply_along_axis(evenFunc, 1, z)
		return div

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
