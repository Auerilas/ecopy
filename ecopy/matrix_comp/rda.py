import numpy as np
from pandas import DataFrame, get_dummies
import matplotlib.pyplot as plt

class rda(object):
	"""
	Docstring for function ecopy.rda
	====================
	Conducts RDA analysis for an site x species matrix Y and
		matrix of predictors X

	Use
	----
	rda(Y, X, scale_y=True, scale_x=False, design_x=False, 
		varNames_y=None, varNames_x=None, rowNames=None,
		pTypes=None)

	Returns an object of class rda

	Parameters
	----------
	Y: A numpy.ndarray or pandas.DataFrame containing species 
		abundances at each site (site x species).
	X:  A pandas.DataFrame or numpy.ndarray containing predictors.
	scale_y: Whether to standardize Y by columns (True) or not (False)
	scale_x: Whether to standardize X by columns (True) or not (False)
	design_x: Whether X has already been made a design matrix (True) or not (False).
	varNames_y: Optional set of variable (i.e. species) names for Y. If None, then
		names are taken from DataFrame columns.
	varNames_x: Optional set of variable (i.e. environmental) names for X. If None, then
		names are taken from DataFrame columns.
	rowNames: Optional set of row (i.e. site) names for Y and X. If None, than 
		names are taken from DataFrame index.
	pTypes: A list denoting whether variables in X are quantitative ('q') or
		factors ('f').

	Attributes
	---------
	spScores: Species scores along each RDA axis
	linSites: Site constraints (linear combinations of predicted values)
	siteScores: Site scores along each RDA axis
	predScores: Scores for each predictor along each RDA axis
	RDA_evals: Eigenvalues for each RDA axis
	corr: Correlation between predictors and each RDA axis
	resid_evals: Eigenvalues for residuals from RDA (based on PCA)
	resid_spScores: Species scores of residuals from RDA (based on PCA)
	resid_siteScores: Site scores of residuals from RDA (based on PCA)
	imp: Summary of importance measures (Std Dev, Proportion Variance) for each
		RDA axis

	Methods
	--------
	summary(): Prints a summary output
	anova(): Conducts an ANOVA to determine significance of the overall model
	triplot(xax=1, yax=2):
		Creates a triplot of species scores, site scores,
			and predictor variable loadings. If predictors are factors,
			they are represented by points. Quantitative predictors
			are represented by arrows.

			xax: which component to put on the x-axis
			yax: which component to put on the y-axis

	Example
	--------
	import ecopy as ep

	dune = ep.load_data('dune')
	dune_env = ep.load_data('dune_env')

	RDA = ep.rda(dune, dune_env[['A1', 'Management']])
	print(RDA.summary())
	print(RDA.anova())
	RDA.triplot()
	"""
	def __init__(self, Y, X, scale_y=True, scale_x=False, design_x=False, varNames_y=None, varNames_x=None, rowNames=None, pTypes=None, sig=False):
		tolerance = 1E-6
		if not isinstance(Y, (DataFrame, np.ndarray)):
			msg = 'Matrix Y must be a pandas.DataFrame or numpy.ndarray'
			raise ValueError(msg)
		if not isinstance(X, (DataFrame, np.ndarray)):
			msg = 'Matrix X must be a pandas.DataFrame or numpy.ndarray'
			raise ValueError(msg)
		if isinstance(X, DataFrame):
			if X.isnull().any().any():
				msg = 'Matrix X contains null values'
				raise ValueError(msg)
		if isinstance(X, np.ndarray):
			if X.dtype=='object':
				msg = 'Matrix X cannot be a numpy.ndarray with object dtype'
				raise ValueError(msg)
			if np.isnan(X).any():
				msg = 'Matrix X contains null values'
				raise ValueError(msg)
		if isinstance(Y, DataFrame):
			if Y.isnull().any().any():
				msg = 'Matrix Y contains null values'
				raise ValueError(msg)
			if (Y.dtypes == 'object').any():
				msg = 'Matrix Y can only contain numeric values'
				raise ValueError(msg)
		if isinstance(Y, np.ndarray):
			if np.isnan(Y).any():
				msg = 'Matrix Y contains null values'
				raise ValueError(msg)
		if varNames_y is None:
			if isinstance(Y, DataFrame):
				varNames_y = Y.columns
			elif isinstance(Y, np.ndarray):
				varNames_y = ['Sp {0}'.format(x) for x in range(1, Y.shape[1]+1)]
		if varNames_x is None:
			if isinstance(X, DataFrame):
				varNames_x = X.columns
			elif isinstance(X, np.ndarray):
				varNames_x = ['Pred {0}'.format(x) for x in range(1, X.shape[1]+1)]
		if rowNames is None:
			if isinstance(Y, DataFrame):
				rowNames = Y.index.values
			elif isinstance(Y, np.ndarray):
				rowNames = ['Site {0}'.format(x) for x in range(1, Y.shape[0]+1)]
		self.y_mat = np.array(Y, dtype='float')
		if design_x:
			self.x_mat = np.array(X, dtype='float')
			if pTypes is None:
				pTypes = ['q']*self.x_mat.shape[1]
		if not design_x and isinstance(X, DataFrame):
			self.x_mat, varNames_x, pTypes = dummyMat(X, scale_x)
		elif not design_x and isinstance(X, np.ndarray):
			self.x_mat = np.array(X, dtype='float')
			if pTypes is None:
				pTypes = ['q']*self.x_mat.shape[1]
			msg = 'Warning: X is a numpy array but not a design matrix. Make sure that matrix X represents the model you wish to analyze. Use patsy.dmatrix if unsure'
			print(msg)
		if not set(pTypes).issubset(['f', 'q']):
			msg = 'pTypes must contain only f or q characters'
			raise ValueError(msg)
		self.y_mat = np.apply_along_axis(lambda x: x-x.mean(), 0, self.y_mat)
		if scale_y:
			self.y_mat = np.apply_along_axis(lambda x: x /x.std(ddof=1), 0, self.y_mat)
		B = np.linalg.pinv(self.x_mat.T.dot(self.x_mat)).dot(self.x_mat.T).dot(self.y_mat)
		yhat = self.x_mat.dot(B)
		yhat_cov = np.cov(yhat, rowvar=0)
		evals, evecs = np.linalg.eig(yhat_cov)
		evals = np.real(evals)
		idx = evals.argsort()[::-1]
		evals = evals[idx]
		evecs = evecs[:,idx]
		RDA_evals = evals[evals>tolerance]
		U = np.real(evecs[:,evals>tolerance])
		F = self.y_mat.dot(U)
		Z = yhat.dot(U)
		C = B.dot(U)
		self.spScores = DataFrame(U.dot(np.diag(RDA_evals**0.5)), index=varNames_y)
		self.linSites = DataFrame(Z.dot(np.diag(RDA_evals**-0.5)), index=rowNames)
		self.siteScores = DataFrame(F.dot(np.diag(RDA_evals**-0.5)), index=rowNames)
		self.predScores = DataFrame(C, index=varNames_x)
		RDA_names = ['RDA Axis {0}'.format(x) for x in range(1, len(RDA_evals)+1)]
		self.spScores.columns=RDA_names
		self.linSites.columns=RDA_names
		self.siteScores.columns=RDA_names
		self.predScores.columns=RDA_names
		self.RDA_evals = RDA_evals
		self.pTypes = pTypes
		self.corr = np.zeros((self.x_mat.shape[1], len(RDA_evals)))
		for i in range(self.x_mat.shape[1]):
			for j in range(len(RDA_evals)):
				self.corr[i,j] = np.corrcoef(self.x_mat[:,i], self.linSites.iloc[:,j])[0,1]
		self.corr = DataFrame(self.corr, index=varNames_x)
		self.corr.columns = RDA_names
		SSY = np.sum(self.y_mat**2)
		SSYhat = np.sum(yhat**2)
		self.R2 = SSYhat / SSY
		self.n = self.y_mat.shape[0]
		self.m = len(B)
		self.R2a = 1. - (1. - self.R2)*((self.n-1.)/(self.n-self.m-1.))
		
		residuals = yhat - self.y_mat
		res_cov = np.cov(residuals, rowvar=0)
		res_evals, res_evecs = np.linalg.eig(res_cov)
		res_evals = np.real(res_evals[res_evals.argsort()[::-1]])
		res_evecs = np.real(res_evecs[:,res_evals.argsort()[::-1]])
		self.resid_evals = res_evals[res_evals > tolerance]
		residU = res_evecs[:,res_evals > tolerance]
		self.resid_spScores = DataFrame(residU.dot(np.diag(self.resid_evals**0.5)), index=varNames_y)
		self.resid_siteScores = DataFrame(residuals.dot(residU).dot(np.diag(self.resid_evals**0.5)), index=rowNames)
		PC_names = ['PC Axis {0}'.format(x) for x in range(1, len(self.resid_evals)+1)]
		self.resid_spScores.columns = PC_names
		self.resid_siteScores.columns = PC_names

		totevals = np.append(self.RDA_evals, self.resid_evals)			 
		sds = np.sqrt(totevals)
		props = totevals / totevals.sum()
		cums = np.cumsum(totevals) / totevals.sum()
		RDA_names.extend(PC_names)
		self.imp = DataFrame(np.vstack((sds, props, cums)), index = ['Std Dev', 'Prop Var', 'Cum Var'])
		self.imp.columns = RDA_names



	def summary(self, n=10):
		print('\nTotal Variance = {0:.3}'.format(np.sum(self.RDA_evals) + np.sum(self.resid_evals)))
		print('Constrained Variance = {0:.3}'.format(np.sum(self.RDA_evals)))
		print('Residual Variance = {0:.3}'.format(np.sum(self.resid_evals)))
		print('R2 = {0:.3}'.format(self.R2))
		print('Adjusted R2 = {0:.3}'.format(self.R2a))

	def anova(self, nperm=999):
		constrained = np.sum(self.RDA_evals)
		resid = np.sum(self.resid_evals)
		Fobs = (constrained/len(self.RDA_evals)) / (resid/len(self.resid_evals))
		Fperm = np.empty(nperm)
		for i in range(nperm):
			idx = np.random.choice(self.n, self.n)
			y_perm = self.y_mat[idx,:]
			B = np.linalg.pinv(self.x_mat.T.dot(self.x_mat)).dot(self.x_mat.T).dot(y_perm)
			yhat_perm = self.x_mat.dot(B)
			cov_perm = np.cov(yhat_perm, rowvar=0)
			evals_perm, evecs_perm = np.linalg.eig(cov_perm)
			evals_perm = np.real(evals_perm)
			evals_perm = evals_perm[evals_perm.argsort()[::-1]]
			evals_perm = evals_perm[evals_perm>1E-6]
			residuals = yhat_perm - self.y_mat
			res_cov = np.cov(residuals, rowvar=0)
			res_evals, res_evecs = np.linalg.eig(res_cov)
			res_evals = np.real(res_evals[res_evals.argsort()[::-1]])
			res_evals[res_evals > 1E-6]
			constrained = np.sum(evals_perm)
			resid = np.sum(res_evals)
			Fperm[i] = (constrained / len(evals_perm)) / (resid/len(res_evals))
		print('Model F-statistic = {0:.3}'.format(Fobs))
		print('p = {0:.4}'.format(np.mean(Fperm > Fobs)))

	def triplot(self, xax=1, yax=2):
		xplot = xax-1
		yplot = yax-1
		f, ax = plt.subplots()
		for i in range(self.spScores.shape[0]):
			ax.plot(self.spScores.iloc[i,xplot], self.spScores.iloc[i,yplot], ms=0)
			ax.text(self.spScores.iloc[i,xplot], self.spScores.iloc[i,yplot], self.spScores.index.values[i], color='r', ha='center', va='center')
		for i in range(self.siteScores.shape[0]):
			ax.plot(self.siteScores.iloc[i,xplot], self.siteScores.iloc[i,yplot], ms=0)
			ax.text(self.siteScores.iloc[i,xplot], self.siteScores.iloc[i,yplot], self.siteScores.index.values[i], color='k', ha='center', va='center')
		for i in range(self.corr.shape[0]):
			if self.pTypes[i] == 'q':
				ax.arrow(0, 0, self.corr.iloc[i,xplot], self.corr.iloc[i,yplot], color='b', head_width=0.1)
				ax.text(self.corr.iloc[i,xplot]*1.2, self.corr.iloc[i,yplot]*1.2, self.corr.index.values[i], color='b', ha='center', va='center')
			elif self.pTypes[i] == 'f':
				ax.text(self.corr.iloc[i,xplot], self.corr.iloc[i,yplot], self.corr.index.values[i], color='b', ha='center', va='bottom')
		ax.set_xlabel('RDA {0}'.format(xax))
		ax.set_ylabel('RDA {0}'.format(yax))
		plt.show()

def dummyMat(matrix, scale):
	rows = matrix.shape[0]
	columns = matrix.shape[1]
	columnType = np.array([""]*columns)
	for col in range(columns):
		id1 = 'q'
		if matrix.iloc[:,col].dtype=='int':
			matrix.iloc[:,col] = matrix.iloc[:,col].apply(float)
		elif matrix.iloc[:,col].dtype=='object':
			id1 = 'f'
		columnType[col] = id1
	modMat = np.array([0]*rows).reshape(rows, 1)
	modNames = ['0']
	pType = ['0']
	for col in range(columns):
		if columnType[col]=='q':
			w = matrix.iloc[:,col]
			w = (w - np.mean(w))
			if scale:
				w = w / np.std(w, ddof=1)
			modMat = np.hstack((modMat, w.reshape(rows, 1)))
			modNames.append(matrix.columns[col])
			pType.append('q')
		elif columnType[col]=='f':
			colName = matrix.columns[col]
			w = get_dummies(matrix.iloc[:,col])
			levels = [colName + ':' + str(x) for x in w.columns]
			modMat = np.hstack((modMat, w))
			modNames.extend(levels)
			pType.extend(['f']*len(levels))
	return modMat[:,1:], modNames[1:], pType[1:]
