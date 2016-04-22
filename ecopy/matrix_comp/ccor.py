import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt

class ccor(object):
	"""
	Docstring for function ecopy.ccor
	====================
	Conducts canonical correlation analysis for two matrices
		Y1 and Y2

	Use
	----
	ccor(Y1, Y2, stand_1=False, stand_2=False, varNames_1=None, varNames_2=None, siteNames=None)

	Returns an object of class ccor

	Parameters
	----------
	Y1: A pandas.DataFrame or numpy.ndarray containing one set of response variables
	Y2: A pandas.DataFrame or numpy.ndarray containing a second set of response variables
	varNames_1: A list of variable names for matrix Y1. If None, inherits from column names of Y1
	varNames_2: A list of variable names for matrix Y2. If None, inherits from column names of Y2
	stand_1: Whether or not to standardize columns of Y1
	stand_2: Whether or not to standardize columns of Y2
	siteNames: A list of site/row names. If None, inherits from index of Y1

	Attributes
	---------
	Scores1: Scores of each row of matrix Y1
	Scores2: Scores of each row of matrix Y2
	loadings1: Variable loadings of Y1
	loadings2: Variable loadings of Y2
	evals: Eigenvalues associated with each axis

	Methods
	--------
	summary(): A summary table for each canonical axis
	biplot(matrix=1, xax=1, yax=1)
		matrix: Which matrix should be plotted. matrix=1 plots the scores and loadings
			of matrix=2, while matrix=2 plots the scores and loadings of matrix 2
		xax: Which canonical axis should on the x-axis
		yax: Which canonical axis should be on the y-axis

	Example
	--------
	import ecopy as ep
	import numpy as np

	Y1 = np.random.normal(size=20*5).reshape(20, 5)
	Y2 = np.random.normal(size=20*3).reshape(20, 3)

	cc = ep.ccor(Y1, Y2)
	cc.summary()
	cc.biplot()
	"""
	def __init__(self, Y1, Y2, varNames_1=None, varNames_2=None, stand_1=False, stand_2=False, siteNames=None):
		if not isinstance(Y1, (DataFrame, np.ndarray)):
			msg = 'Matrix Y1 must be a pandas.DataFrame or numpy.ndarray'
			raise ValueError(msg)
		if not isinstance(Y2, (DataFrame, np.ndarray)):
			msg = 'Matrix Y2 must be a pandas.DataFrame or numpy.ndarray'
			raise ValueError(msg)
		if isinstance(Y2, DataFrame):
			if Y2.isnull().any().any():
				msg = 'Matrix Y2 contains null values'
				raise ValueError(msg)
		if isinstance(Y2, np.ndarray):
			if Y2.dtype=='object':
				msg = 'Matrix Y2 cannot be a numpy.ndarray with object dtype'
				raise ValueError(msg)
			if np.isnan(Y2).any():
				msg = 'Matrix Y2 contains null values'
				raise ValueError(msg)
		if isinstance(Y1, DataFrame):
			if Y1.isnull().any().any():
				msg = 'Matrix Y1 contains null values'
				raise ValueError(msg)
			if (Y1.dtypes == 'object').any():
				msg = 'Matrix Y1 can only contain numeric values'
				raise ValueError(msg)
		if isinstance(Y1, np.ndarray):
			if np.isnan(Y1).any():
				msg = 'Matrix Y1 contains null values'
				raise ValueError(msg)
		if varNames_1 is None:
			if isinstance(Y1, DataFrame):
				varNames_1 = Y1.columns
			elif isinstance(Y1, np.ndarray):
				varNames_1 = ['Y1 {0}'.format(x) for x in range(1, Y1.shape[1]+1)]
		if varNames_2 is None:
			if isinstance(Y2, DataFrame):
				varNames_2 = Y2.columns
			elif isinstance(Y2, np.ndarray):
				varNames_2 = ['Y2 {0}'.format(x) for x in range(1, Y2.shape[1]+1)]
		if siteNames is None:
			if isinstance(Y1, DataFrame):
				siteNames = Y1.index.values
			elif isinstance(Y1, np.ndarray):
				siteNames = ['Site {0}'.format(x) for x in range(1, Y1.shape[0]+1)]
		if Y1.shape[0] != Y2.shape[0]:
			msg = 'Matrices must have same number of rows'
			raise ValueError(msg)
		Y1 = np.array(Y1)
		Y2 = np.array(Y2)
		if stand_1:
			Y1 = (Y1 - Y1.mean(axis=0)) / Y1.std(axis=0)
		if stand_2:
			Y2 = (Y2 - Y2.mean(axis=0)) / Y2.std(axis=0)
		df = float(Y1.shape[0] - 1)
		D1 = Y1 - Y1.mean(axis=0)
		D2 = Y2 - Y2.mean(axis=0)
		S1 = D1.T.dot(D1) * 1./df
		S2 = D2.T.dot(D2) * 1./df
		S12 = D1.T.dot(D2) * 1./df
		Chol1 = np.linalg.pinv(np.linalg.cholesky(S1).T)
		Chol2 = np.linalg.pinv(np.linalg.cholesky(S2).T)
		K = Chol1.T.dot(S12).dot(Chol2)
		V, W, U = np.linalg.svd(K)
		U = U.T
		CoefY1 = Chol1.dot(V)
		CoefY2 = Chol2.dot(U)
		self.Scores1 = DataFrame(Y1.dot(CoefY1), index=siteNames)
		self.Scores2 = DataFrame(Y2.dot(CoefY2), index=siteNames)
		self.loadings1 = np.corrcoef(Y1, self.Scores1, rowvar=0)[:5, 5:]
		self.loadings2 = np.corrcoef(Y2, self.Scores2, rowvar=0)[:3, 3:]
		axes1 = ['CA Axis {0}'.format(x) for x in range(1, self.loadings1.shape[1]+1)]
		self.loadings1 = DataFrame(self.loadings1, index=varNames_1, columns=axes1)
		axes2 = ['CA Axis {0}'.format(x) for x in range(1, self.loadings2.shape[1]+1)]
		self.loadings2 = DataFrame(self.loadings2, index=varNames_2, columns=axes2)
		self.evals = W
		

	def summary(self):
		print('Constrained variance = {0:.3}'.format(np.sum(self.evals)))
		print('Constrained variance explained be each axis')
		print([str(i) for i in np.round(self.evals, 3)])
		print('Proportion constrained variance')
		print([str(i) for i in np.round(self.evals/self.evals.sum(), 3)])

	def biplot(self, matrix=1, xax=1, yax=2):
		xplot = xax-1
		yplot = yax-1
		scores = self.Scores1
		loadings = self.loadings1
		if matrix==2:
			scores = self.Scores2
			loadings = self.loadings2
		varNames = loadings.index.values
		siteNames = scores.index.values
		f, ax = plt.subplots()
		for i in range(scores.shape[0]):
			ax.plot(scores.iloc[i,xplot], scores.iloc[i,yplot], ms=0)
			ax.text(scores.iloc[i,xplot], scores.iloc[i,yplot], siteNames[i], color='r', ha='center', va='center')
		for i in range(loadings.shape[0]):
			ax.arrow(0, 0, loadings.iloc[i, xplot], loadings.iloc[i, yplot], color='b', head_width=0.1)
			ax.text(loadings.iloc[i,xplot]*1.2, loadings.iloc[i,yplot]*1.2, varNames[i], color='b', ha='center', va='center')
		xmins = (scores.iloc[:,xplot].min(), loadings.iloc[:,xplot].min()*1.2)
		xmax = (scores.iloc[:,xplot].max(), loadings.iloc[:,xplot].max()*1.2)
		ymins = (scores.iloc[:,yplot].min(), loadings.iloc[:,yplot].min()*1.2)
		ymax = (scores.iloc[:,yplot].max(), loadings.iloc[:,yplot].max()*1.2)
		ax.set_xlabel('CCor {0}'.format(xax))
		ax.set_ylabel('CCor {0}'.format(yax))
		ax.set_xlim([min(xmins), max(xmax)])
		ax.set_ylim([min(ymins), max(ymax)])
		plt.show()
		




