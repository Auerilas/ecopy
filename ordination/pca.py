import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as py
import warnings

class pca:
	'''
	Docstring for function ecopy.pca
	====================
	Conducts principle components analysis (PCA). User
		can supply either a correlation matrix, covariance matrix,
		or dataframe containing only the variables for ordination

	Use
	----
	pca(self, x, scaled = True, varNames = None)

	Returns an object of class pca

	Parameters
	----------
	x:  Data for ordination. Should either be a pandas DataFrame or numpy array
	scale: Boolean for whether the data should be standardized prior to pca
	varNames: If x is numpy array, pass a list of variable names corresponding to the columns

	Attributes
	---------
	evals: array of eigenvalues
	evecs: array of eigenvectors
	scores: array of observation scores along each component

	Methods
	--------
	summary(): print a summary table
	biplot(scale = 0, obsNames = False): create a biplot of the first two components with a given scaling
		scales each score by 1/lamba**scale (shrinking proportional to sd of component)
		scales each vector by lambda**scale (expanding vector proportional to sd of component)
		scale should be between 0 and 1
		names tells whether to plot the names of each observation (False by default)

	Example
	--------
	import pandas.rpy.common as com
	from ecopy import pca
	USArrests = com.load_data('USArrests')
	prcomp = pca(USArrests, scaled = True)
	prcomp.summary()
	prcomp.biplot(scale = 0)
	prcomp.biplot(scale = 1, obsNames = True)
	'''
	def __init__(self, x, scale = True, varNames = None):
		# if the data is not a dataframe or array, raise error
		if not isinstance(x, (DataFrame, np.ndarray)):
			msg = 'Data must either be DataFrame or array'
			raise ValueError(msg)
		#  if x is a DataFrame
		if isinstance(x, DataFrame):
			# check NAs
			if x.isnull().any().any():
				msg = 'DataFrame contains null values'
				raise ValueError(msg)
			# check for non-numeric
			if (x.dtypes == 'object').any():
				msg = 'DataFrame can only contain numeric values'
				raise ValueError(msg)
			# convert to a numpy array
			y = np.array(x)
		# if x is array, simple re-assign	
		if isinstance(x, np.ndarray):
			if np.isnan(x).any():
				msg = 'Array contains null values'
				raise ValueError(msg)
			y = x
		# check that scale is boolean
		if not isinstance(scale, bool):
			msg = "scale argument must be boolean"
			raise ValueError(msg)
		# if scale, standardize data
		if scale:
			y = np.apply_along_axis(standardize, 0, y)
		# get covariance matrix
		covMat = np.cov(y, rowvar = 0)
		self.evals, self.evecs = np.linalg.eig(covMat)
		idx = self.evals.argsort()[::-1]
		self.evals = self.evals[idx]
		self.evecs = self.evecs[:,idx]
		self.scores = y.dot(self.evecs)
		# save the output into dataframes
		names = ['PC' + str(i) for i in range(1, self.evecs.shape[1] + 1)]
		self.sd = np.sqrt(self.evals)
		self.props = self.evals/np.sum(self.evals)
		self.imp = DataFrame(np.vstack((self.sd, self.props)), index = ['Std Dev', 'Cum Prop'])
		self.imp.columns = names
		if isinstance(x, DataFrame):
			self.rot = DataFrame(self.evecs, index = x.columns)
			self.rot.columns = names
		else:
			self.rot = DataFrame(self.evecs, index = varNames)
			self.rot.columns = names
		if isinstance(x, DataFrame):
			self.labs =  np.array(x.index)
		else:
			self.labs = range(x.shape[0])

	def summary(self):
		print 
		print 'Component Importance:'
		print self.imp
		print
		print 'Rotation:'
		print self.rot

	def biplot(self, scale = 0, obsNames = False):
		if scale > 1:
			msg = 'Warning: Scale should be in the range [0,1]'
			warnings.warn(msg)
		# shrink points proportional to the inverse standard deviation of each component	
		ScorePlot = self.scores.dot(np.diag(1./np.sqrt(self.evals)**scale))
		# expand vectors proportion to standard deviation of each component
		VecPlot = self.evecs.dot(np.diag(np.sqrt(self.evals)**scale))
		py.scatter(ScorePlot[:,0], ScorePlot[:,1])
		if obsNames:
			for i in range(ScorePlot.shape[0]):
				py.text(ScorePlot[i,0], ScorePlot[i,1], self.labs[i], ha = 'center', va = 'center')
		for i in range(VecPlot.shape[0]):
			py.arrow(0, 0, VecPlot[i,0], VecPlot[i,1], color = 'red')
			py.text(VecPlot[i, 0]*1.1, VecPlot[i,1]*1.1, self.rot.index[i], color = 'red', ha = 'center', va = 'center')
		py.xlabel('PC 1')
		py.ylabel('PC 2')
		py.show()



def standardize(a):
	return (a - np.mean(a))/np.std(a, ddof = 1)