import numpy as np
from pandas import DataFrame 

def transform(x, method='wisconsin', axis=1, breakNA=True):
	'''
	Docstring for function ecopy.distance
	========================
	Applies a transformation to a given matrix

	Use
	----
	transform(x, method='wisconsin', axis=1, breakNA=True)

	Returns either a numpy.ndarray or pandas.DataFrame, depending
		on input

	Parameters
	----------
	x: pandas dataframe or numpy array with observations as rows
		and descriptors as columns
	method: particular transformation
		total: divides by the column or row sum
		max: divides by the column or row max
		normalize: makes the sum of squares = 1 along columns or rows (chord distance)
		range: converts to a range of [0, 1]
		standardize: calculates z-score along columns or rows
		hellinger: square-root after transformation by total
		log: return ln(x + 1) which automatically returns 0 if x = 0
		pa: convert to binary presence/absence
		wisconsin: the standard transformation, first by column maxima then by 
			row totals (default)
	axis: which axis should undergo transformation.
		axis = 0 applies down columns
		axis = 1 applies across rows (default)
	breakNA: should the process halt if the matrix contains any NAs?
		if False, then NA's are ignored during transformations
	
	Example
	--------
	import ecopy as ep
	varespec = ep.load_data('varespec')

	# divide each element by row total
	ep.transform(varespec, method='total', axis=1)
	'''
	if not isinstance(breakNA, bool):
		msg = 'breakNA must be boolean'
		raise ValueError(msg)
	if not isinstance(x, (DataFrame, np.ndarray)):
		msg = 'x must be either numpy array or dataframe'
		raise ValueError(msg)
	if axis not in [0, 1]:
		msg = 'Axis argument must be either 0 or 1'
		raise ValueError(msg)
	if method not in ['total', 'max', 'normalize', 'range', 'standardize', 'hellinger', 'log', 'logp1', 'pa', 'wisconsin']:
		msg = '{0} not an accepted method'.format(method)
		raise ValueError(msg)
	if isinstance(x, DataFrame):
		if breakNA:
			if x.isnull().any().any():
				msg = 'DataFrame contains null values'
				raise ValueError(msg)
		if (x<0).any().any():
			msg = 'DataFrame contains negative values'
			raise ValueError(msg)
		z = x.copy()
		if method=='total':
			data = z.apply(totalTrans, axis=axis)
			return data
		if method=='max':
			data = z.apply(maxTrans, axis=axis)
			return data
		if method=='normalize':
			data = z.apply(normTrans, axis=axis)
			return data
		if method=='range':
			data = z.apply(rangeTrans, axis=axis)
			return data
		if method=='standardize':
			data = z.apply(standTrans, axis=axis)
			return data
		if method=='pa':
			z[z>0] = 1
			return z
		if method=='hellinger':
			data = z.apply(totalTrans, axis=axis)
			return np.sqrt(data)
		if method=='log':
			if ((z > 0) & (z < 1)).any().any():
				msg = 'Log of values between 0 and 1 will return negative numbers\nwhich cannot be used in subsequent distance calculations'
				raise ValueError(msg)
			data = z.applymap(lambda y: np.log(y+1))
			return data
		if method=='logp1':
			if ((z > 0) & (z < 1)).any().any():
				msg = 'Log of values between 0 and 1 will return negative numbers\nwhich cannot be used in subsequent distance calculations'
				raise ValueError(msg)
			data = z.applymap(lambda y: np.log(y) + 1 if y>0 else 0)
			return data
		if method=='wisconsin':
			data = z.apply(maxTrans, axis=0)
			data = data.apply(totalTrans, axis=1)
			return data
	if isinstance(x, np.ndarray):
		if breakNA:
			if np.isnan(np.sum(x)):
				msg = 'Array contains null values'
				raise ValueError(msg)
		if (x < 0).any():
			msg = 'Array contains negative values'
			raise ValueError(msg)
		z = x.copy()
		if method=='total':
			data = np.apply_along_axis(totalTrans, axis, z)
			return data
		if method=='max':
			data = np.apply_along_axis(maxTrans, axis, z)
			return data
		if method=='normalize':
			data = np.apply_along_axis(normTrans, axis, z)
			return data
		if method=='range':
			data = np.apply_along_axis(rangeTrans, axis, z)
			return data
		if method=='standardize':
			data = np.apply_along_axis(standTrans, axis, z)
			return data
		if method=='pa':
			z[z>0] = 1
			return z
		if method=='hellinger':
			data = np.apply_along_axis(totalTrans, axis, z)
			return np.sqrt(data)
		if method=='log':
			if ((z > 0) & (z < 1)).any():
				msg = 'Log of values between 0 and 1 will return negative numbers\nwhich cannot be used in subsequent distance calculations'
				raise ValueError(msg)
			data = np.log(z + 1)
			return data
		if method=='logp1':
			if ((z > 0) & (z < 1)).any():
				msg = 'Log of values between 0 and 1 will return negative numbers\nwhich cannot be used in subsequent distance calculations'
				raise ValueError(msg)
			data = z.astype('float')
			data[np.greater(data,0)] = np.log(data[np.greater(data,0)]) + 1
			return data
		if method=='wisconsin':
			data = np.apply_along_axis(maxTrans, 0, z)
			data = np.apply_along_axis(totalTrans, 1, z)
			return data

def totalTrans(y):
	return y/np.nansum(y)

def maxTrans(y):
	return y/np.nanmax(y)

def normTrans(y):
	denom = np.nansqrt(np.nansum(y**2))
	return y/denom

def rangeTrans(y):
	return (y - np.nanmin(y))/(np.nanmax(y) - np.nanmin(y))

def standTrans(y):
	return (y - np.nanmean(y))/np.nanstd(y, ddof=1)
