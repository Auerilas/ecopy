import numpy as np
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
from .isoFunc import _isotonic_regression

class isotonic:
	'''
	Docstring for function ecopy.isotonic
	=============================

	Use
	----
	isotonic(y, x=None, w=None, direction='increasing')

	Returns an object of class isotonic

	Parameters
	----------------
	y: A list, pandas.Series, pandas.DataFrame, or numpy.ndarray
		giving the response variable
	x:  A list, pandas.Series, pandas.DataFrame, or numpy.ndarray
		giving the predictor variable (optional)
	w:  A list, pandas.Series, pandas.DataFrame, or numpy.ndarray
		giving observation weights (optional)
	direction: Specify whether the function should be
		'increasing' or 'decreasing'

	Attributes
	---------
	prediction: An np.nadarray of the predicted values of each observation
	predictor: The original predictor variable (x)
	obs: The observations (y)
	
	Methods
	--------
	plot(): A plot of the original observations and the fitted values
		against the predictor values (if specified)
	
	Example
	--------
	import numpy as np
	import ecopy as ep
	age = np.array([8,  8,  8, 10, 10, 10, 12, 12, 12, 14, 14])
	pit = np.array([21.0, 23.5, 23.0, 24.0, 21.0, 25.0, 21.5, 22.0, 19.0, 23.5, 25.0])
	solution = ep.isotonic(pit, age)
	solution.plot()
	'''
	def __init__(self, y, x=None, w=None, direction='increasing'):
		if w is None:
			w = np.ones(len(y))
		if x is None:
			x = np.arange(len(y))
		if not isinstance(y, (DataFrame, Series, np.ndarray, list)):
			msg = 'Response variable (y) must be a pandas.DataFrame, pandas.Series, numpy.ndarray, or list'
			raise ValueError(msg)
		if not isinstance(x, (DataFrame, Series, np.ndarray, list)):
			msg = 'Predictor variable (x) must be a pandas.DataFrame, pandas.Series, numpy.ndarray, or list'
			raise ValueError(msg)
		if not isinstance(w, (DataFrame, Series, np.ndarray, list)):
			msg = 'Weights vector (w) must be a pandas.DataFrame, pandas.Series, numpy.ndarray, or list'
			raise ValueError(msg)
		if direction not in ['increasing', 'decreasing']:
			msg = "direction must be either 'increasing' or 'decreasing'"
			raise ValueError(msg)
		x = np.array(x)
		y = np.array(y)
		w = np.array(w)
		o = x.argsort() # sorting
		o2 = np.zeros(len(x), dtype='int')
		for i in range(len(y)):
			o2[o[i]] = i # get the original order (i.e. return the second sorted value to its original position)
		self.predictor =  x
		y = y[o]
		w = w[o]
		if direction is 'decreasing':
			y = y[::-1]	# flip y around and use the same algorithm
		self.prediction = _isotonic_regression(y, w, np.ones(len(y)))
		self.prediction = self.prediction[o2] # put into original order
		self.obs = y[o2]

	def plot(self):
		o = self.predictor.argsort()
		xplt = self.predictor[o]
		yplt = self.obs[o]
		predplt = self.prediction[o]
		f, ax = plt.subplots()
		ax.scatter(xplt, yplt, c='r', s=50, label='Original Data')
		ax.step(xplt, predplt, label='Prediction')
		ax.set_ylabel("Response")
		ax.set_xlabel("Predictor")
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.legend(loc='best')
		plt.show()