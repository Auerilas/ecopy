import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as py

class pca(object):
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

	Methods
	--------
	summary_imp(): print a summary table of principle axes
	summary_rot(): print a summary table of axes rotations
	summary_desc():cumulative variance explained by each principle axis for each
		descriptor. 
	biplot(xax=1, yax=2, type = "distance", obsNames = False):
		create a biplot of the first two components with a given scaling
			xax: which component to put on the x-axis
			yax: which component to put on the y-axis
			type: denotes whether a distance biplot or correlation
				biplot should be used
			obsNames: tells whether to plot the names of each observation (False by default)

	Example
	--------
	import ecopy as ep
	USArrests = ep.load_data('USArrests')
	prcomp = ep.pca(USArrests.iloc[:,1:], scale = True)
	prcomp.summary_imp()
	prcomp.summary_rot()
	prcomp.summary_desc()
	prcomp.biplot(type = 'distance')
	prcomp.biplot(type = 'correlation', obsNames = True)
	'''
	def __init__(self, x, scale = True, varNames = None):
		if not isinstance(x, (DataFrame, np.ndarray)):
			msg = 'Data must either be pandas.DataFrame or numpy.ndarray'
			raise ValueError(msg)
		if isinstance(x, DataFrame):
			if x.isnull().any().any():
				msg = 'DataFrame contains null values'
				raise ValueError(msg)
			if (x.dtypes == 'object').any():
				msg = 'DataFrame can only contain numeric values'
				raise ValueError(msg)
			y = np.array(x)
		if isinstance(x, np.ndarray):
			if np.isnan(x).any():
				msg = 'Array contains null values'
				raise ValueError(msg)
			y = x
		if not isinstance(scale, bool):
			msg = "scale argument must be boolean"
			raise ValueError(msg)
		if scale:
			y = np.apply_along_axis(standardize, 0, y)
		else:
			y = np.apply_along_axis(lambda z: z - np.mean(z), 0, y)
		self.evals, self.evecs, self.scores = eig_decomp(y)
		if isinstance(x, DataFrame):
			self.varNames = x.columns
		else:
			self.varNames = varNames
		if isinstance(x, DataFrame):
			self.labs =  np.array(x.index)
		else:
			self.labs = range(x.shape[0])



	def summary_imp(self):
		sd = np.sqrt(self.evals)
		props = self.evals/np.sum(self.evals)
		cums = np.cumsum(self.evals)/np.sum(self.evals)
		names = ['PC' + str(i) for i in range(1, self.evecs.shape[1]+1)]
		imp = DataFrame(np.vstack((sd, props, cums)), index = ['Std Dev', 'Prop Var', 'Cum Var'])
		imp.columns = names
		return imp

	def summary_rot(self):
		names = ['PC' + str(i) for i in range(1, self.evecs.shape[1]+1)]
		rot = DataFrame(self.evecs, index=self.varNames)
		rot.columns = names
		return rot

	def summary_desc(self):
		names = ['PC' + str(i) for i in range(1, self.evecs.shape[1]+1)]
		U2 = self.evecs.dot(np.diag(self.evals**0.5))
		desCums = np.apply_along_axis(lambda x: np.cumsum(x**2) / np.sum(x**2), 1, U2)
		desc = DataFrame(desCums, index=self.varNames)
		desc.columns = names
		return desc

	def biplot(self, xax=1, yax=2, type='distance', obsNames = False):
		if type not in ['distance', 'correlation']:
			msg = 'type argument must be either distance or correlation'
			raise ValueError(msg)
		if type=='distance':
			ScorePlot = self.scores
			VecPlot = self.evecs
		if type=='correlation':
			VecPlot = self.evecs.dot(np.diag(self.evals**0.5))		
			ScorePlot = self.scores.dot(np.diag(self.evals**-0.5))
		if self.varNames is None:
			self.varNames = range(1, self.evecs.shape[0]+1)
		f, ax = py.subplots()
		ax.axvline(0, ls='solid', c='k')
		ax.axhline(0, ls='solid', c='k')
		if obsNames:
			ax.scatter(ScorePlot[:,xax-1], ScorePlot[:,yax-1], s=0)
			for i in range(ScorePlot.shape[0]):
				py.text(ScorePlot[i,xax-1], ScorePlot[i,yax-1], self.labs[i], ha = 'center', va = 'center')
		else:
			ax.scatter(ScorePlot[:,xax-1], ScorePlot[:,yax-1])
		for i in range(VecPlot.shape[0]):
			ax.arrow(0, 0, VecPlot[i,xax-1], VecPlot[i,yax-1], color = 'red', head_width=.05)
			ax.text(VecPlot[i, xax-1]*1.2, VecPlot[i,yax-1]*1.2, self.varNames[i], color = 'red', ha = 'center', va = 'center')
		xmax = max(np.amax(ScorePlot[:,xax-1]), np.amax(VecPlot[:,xax-1]))
		xmin = min(np.min(ScorePlot[:,xax-1]), np.min(VecPlot[:,xax-1]))
		ymax = max(np.amax(ScorePlot[:,yax-1]), np.amax(VecPlot[:,yax-1]))
		ymin = min(np.amin(ScorePlot[:,yax-1]), np.amin(VecPlot[:,yax-1]))
		ax.set_xlim([xmin + 0.15*xmin, xmax+0.15*xmax])
		ax.set_ylim([ymin + 0.15*ymin, ymax+0.15*ymax])
		ax.set_xlabel('PC {!s}'.format(xax))
		ax.set_ylabel('PC {!s}'.format(yax))
		py.show()



def standardize(a):
	return (a - np.mean(a))/np.std(a, ddof = 1)

def eig_decomp(y):
	V, W, U = np.linalg.svd(y)
	N = y.shape[0]
	evecs = U.T
	evals = W**2/(N-1)
	scores = y.dot(evecs)
	return np.real(evals), np.real(evecs), scores