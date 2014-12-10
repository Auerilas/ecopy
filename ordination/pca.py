import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as py

class pca:
	'''
	Docstring for function ecopy.pca
	====================
	Conducts principle components analysis (PCA). User
		supplies a matrix containing only the variables for ordination

	Use
	----
	pca(x, scale = True, varNames = None)

	Returns an object of class pca

	Parameters
	----------
	x:  Data for ordination. Should either be a pandas DataFrame or numpy.ndarray.
		Observations as rows and descriptors as columns
	scale: Boolean for whether the data should be standardized prior to pca
	varNames: If x is numpy array, pass a list of variable names corresponding to the columns

	Attributes
	---------
	evals: array of eigenvalues
	evecs: array of eigenvectors
	scores: array of observation scores along each component
	correlation: correlation of each descriptor with each component
	cumdesc: cumulative variance explained by each principle axis for each
		descriptor. 

	Methods
	--------
	summary_imp(): print a summary table of principle axes
	summary_rot(): print a summary table of axes rotations
	biplot(type = "distance", obsNames = False):
		create a biplot of the first two components with a given scaling
			type: denotes whether a distance biplot or correlation
				biplot should be used
			obsNames: tells whether to plot the names of each observation (False by default)

	Example
	--------
	import pandas.rpy.common as com
	from ecopy import pca
	USArrests = com.load_data('USArrests')
	prcomp = pca(USArrests, scale = True)
	prcomp.summary()
	prcomp.correlation
	prcomp.biplot(scale = 1)
	prcomp.biplot(scale = 0.5, obsNames = True)
	'''
	def __init__(self, x, scale = True, varNames = None):
		# if the data is not a dataframe or array, raise error
		if not isinstance(x, (DataFrame, np.ndarray)):
			msg = 'Data must either be pandas.DataFrame or nump.ndarray'
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
		else:
			y = np.apply_along_axis(lambda z: z - np.mean(z), 0, y)
		# get covariance matrix
		covMat = np.cov(y, rowvar = 0)
		self.evals, self.evecs = np.linalg.eig(covMat)
		idx = self.evals.argsort()[::-1]
		self.evals = self.evals[idx]
		self.evecs = self.evecs[:,idx]
		self.scores = y.dot(self.evecs)
		# save the output into dataframes
		names = ['PC' + str(i) for i in range(1, self.evecs.shape[1] + 1)]
		sd = np.sqrt(self.evals)
		props = self.evals/np.sum(self.evals)
		cums = np.cumsum(self.evals) / np.sum(self.evals)
		self.imp = DataFrame(np.vstack((sd, props, cums)), index = ['Std Dev', 'Proportion', 'Cum Prop'])
		self.imp.columns = names
		self.correlation = np.zeros([y.shape[1], y.shape[1]])
		U2 = self.evecs.dot(np.diag(self.evals**0.5))
		desCums = np.apply_along_axis(lambda x: np.cumsum(x**2) / np.sum(x**2), 1, U2)
		for i in xrange(y.shape[1]):
			for j in xrange(y.shape[1]):
				self.correlation[i,j] = self.evecs[i,j]*np.sqrt(self.evals[i])/ np.sqrt(covMat[j,j])
		if isinstance(x, DataFrame):
			self.rot = DataFrame(self.evecs, index = x.columns)
			self.rot.columns = names
			self.correlation = DataFrame(self.correlation, index=x.columns)
			self.correlation.columns = names
			self.cumdesc = DataFrame(desCums, index=x.columns)
			self.cumdesc.columns = names
		else:
			self.rot = DataFrame(self.evecs, index = varNames)
			self.rot.columns = names
			self.correlation = DataFrame(self.correlation, index=varNames)
			self.correlation.columns = names
			self.cumdesc = DataFrame(desCums, index=varNames)
			self.cumdesc.columns = names
		if isinstance(x, DataFrame):
			self.labs =  np.array(x.index)
		else:
			self.labs = range(x.shape[0])

	def summary_imp(self):
		return self.imp

	def summary_rot(self):
		return self.rot

	def biplot(self, type='distance', obsNames = False):
		if type not in ['distance', 'correlation']:
			msg = 'type argument must be either distance or correlation'
			raise ValueError(msg)
		if type=='distance':
			ScorePlot = self.scores
			VecPlot = self.evecs
		if type=='correlation':
			VecPlot = self.evecs.dot(np.diag(self.evals**0.5))		
			ScorePlot = self.scores.dot(np.diag(self.evals**-0.5))
		f, ax = py.subplots()
		if obsNames:
			ax.scatter(ScorePlot[:,0], ScorePlot[:,1], s=0)
			for i in range(ScorePlot.shape[0]):
				py.text(ScorePlot[i,0], ScorePlot[i,1], self.labs[i], ha = 'center', va = 'center')
		else:
			ax.scatter(ScorePlot[:,0], ScorePlot[:,1])
		for i in range(VecPlot.shape[0]):
			ax.arrow(0, 0, VecPlot[i,0], VecPlot[i,1], color = 'red', head_width=.05)
			ax.text(VecPlot[i, 0]*1.1, VecPlot[i,1]*1.1, self.rot.index[i], color = 'red', ha = 'center', va = 'center')
		xmax = max(np.amax(ScorePlot[:,0]), np.amax(VecPlot[:,0]))
		xmin = min(np.min(ScorePlot[:,0]), np.min(VecPlot[:,0]))
		ymax = max(np.amax(ScorePlot[:,1]), np.amax(VecPlot[:,1]))
		ymin = min(np.amin(ScorePlot[:,1]), np.amin(VecPlot[:,1]))
		pMin = min(xmin, ymin)
		pMax = max(xmax, ymax)
		ax.set_xlim([pMin - 0.15*pMin, pMax+0.15*pMax])
		ax.set_ylim([pMin - 0.15*pMin, pMax+0.15*pMax])
		ax.set_xlabel('PC 1')
		ax.set_ylabel('PC 2')
		py.show()



def standardize(a):
	return (a - np.mean(a))/np.std(a, ddof = 1)