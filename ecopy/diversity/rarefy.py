import numpy as np
from pandas import DataFrame
from scipy.misc import comb
import matplotlib.pyplot as plt

def rarefy(x, method='rarefy', size = None, breakNA=True):
	'''
	Docstring for function ecopy.rarefy
	========================
	Various rarefaction techniques for a site x species matrix.
		All indices computed along rows (axis = 1)

	Use
	----
	rarefy(x, method='rarefy', size=None, breakNA=True)

	Parameters
	----------
	x:  numpy array or pandas dataframe with observations as rows
		and descriptors as columns
	method: a method used for rarefaction
		rarefy: Calculates estimated richness rarified to a given sample 
			size (see size parameter).
			sum(1 - nCr(N-Ni, size) / (nCr(N, size)))
		rareCurve: Draws a rarefaction curve for each site (row). 
			Rarefaction curves use the following functoin
			Sn - sum(1 - nCr(N-Ni, i))/nCr(N, i)
			
	size: the sample size used in rarefaction. Can be left empty,
		in which case size is the minimum of row sums (number of 
		individuals from the sparsest site). Can be a single number, which applies
		the same size to all rows. Can be a numpy array that contains 
		different sizes for each site.
	breakNA: should the process halt if the matrix contains any NAs?
		if False, then NA's undergo pairwise deletion during distance calculation,
		such that when calculating the distance between two rows, if any
		species is missing from a row, that species is removed from both rows

	Example
	--------
	import ecopy as ep
	BCI = ep.load_data('BCI')
	
	# calculate rarefied species richness
	rareRich = ep.rarefy(BCI, 'rarefy')

	# draw rarefaction curves
	ep.rarefy(BCI, 'rarecurve')
	'''
	listofmethods = ['rarefy', 'rarecurve']
	if not isinstance(breakNA, bool):
		msg = 'removaNA argument must be boolean'
		raise ValueError(msg)
	if method not in listofmethods:
		msg = 'method argument {0!s} is not an accepted rarefaction method'.format(method)
		raise ValueError(msg)
	if not isinstance(x, (DataFrame, np.ndarray)):
		msg = 'x argument must be a numpy array or pandas dataframe'
		raise ValueError(msg)
	if size is not None:
		if not isinstance(size, (int, float, np.ndarray)):
			msg = 'size must be integer, float, or numpy array'
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
		if method=='rarefy':
			if size is None:
				sums = x.apply(sum, axis=1)
				size = np.min(sums)
				rich = x.apply(rare, axis=1, args=(size,))
				return rich
			else:
				if isinstance(size, (int, float)):
					 rich = x.apply(rare, axis=1, args=(size,))
					 return rich
				else:
					if len(size) != len(x):
						msg = 'length of size does not match number of rows'
						raise ValueError(msg)
					z = x.copy()
					z['size'] = size
					rich = z.apply(rare_wrapper, axis=1)
					return rich
		if method=='rarecurve':
			z = x.copy()
			z.reset_index(inplace=True)
			z.apply(rCurve, axis=1)
			plt.xlabel('Number of Individuals')
			plt.ylabel('Number of Species')
			plt.show()
	if isinstance(x, np.ndarray):
		if breakNA:
			if np.isnan(np.sum(x)):
				msg = 'Array contains null values'
				raise ValueError(msg)
		if (x < 0).any():
			msg = 'Array contains negative values'
			raise ValueError(msg)
		if method=='rarefy':
			if size is None:
				sums = np.apply_along_axis(np.nansum, 1, x)
				size = np.min(sums)
				rich = np.apply_along_axis(rare, 1, x, size)
				return rich
			else:
				if isinstance(size, (int, float)):
					 rich = np.apply_along_axis(rare, 1, x, size)
					 return rich
				else:
					if len(size) != x.shape[0]:
						msg = 'length of size does not match number of rows'
						raise ValueError(msg)
					N = np.nansum(x, axis=1)
					diff = (N[:,np.newaxis] - x).T
					return np.sum(1 - comb(diff, size)/comb(N, size), axis=0)
		if method=='rarecurve':
			z = DataFrame(x)
			z.reset_index(inplace=True)
			z.apply(rCurve, axis=1)
			plt.xlabel('Number of Individuals')
			plt.ylabel('Number of Species')
			plt.show()

def rare(y, size):
	notabs = ~np.isnan(y)
	t = y[notabs]
	N = np.sum(t)
	diff = N - t
	rare_calc = np.sum(1 - comb(diff, size)/comb(N, size))
	return rare_calc

def rare_wrapper(data):
	s2 = data['size']
	x2 = data.drop('size')
	return rare(x2, s2)

def rareCurve_Func(i, Sn, n, x):
	sBar = Sn -  np.sum(comb(n-x, i))/comb(n, i) 
	return sBar

def rCurve(x):
	ix = x['index']
	z = x.drop('index').astype('float')
	notabs = ~np.isnan(z)
	y = z[notabs]
	n = np.sum(y)
	Sn = len(z)
	iPred = np.linspace(0, n, 1000)
	yhat = [rareCurve_Func(i, Sn, n, y) for i in iPred]
	plt.plot(iPred, yhat)
	plt.text(iPred[-1], yhat[-1], str(ix), ha='left', va='center')
