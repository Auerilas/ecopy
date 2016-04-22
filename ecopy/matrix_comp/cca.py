import numpy as np
from pandas import DataFrame, get_dummies
import matplotlib.pyplot as plt
from ecopy import ca

class cca(object):
	"""
	Docstring for function ecopy.cca
	====================
	Conducts canonical correspondance analysis (CA) for a given
		response matrix (Y) and predictor variables (X)

	Use
	----
	cca(Y, X, varNames_y=None, varNames_x=None, rowNames=None, scaling=1)

	Returns an object of class cca

	Parameters
	----------
	Y: A pandas.DataFrame or numpy.ndarray containing species abundance data (site x species)
	X: A pandas.DataFrame or numpy.ndarray containing predictor variables for constrained ordination (site x variable).
	varNames_y: A list of variable names for matrix Y. If None, inherits from column names of Y.
	varNames_x: A list of variable names for matrix X. If None, inherits from column names of X.
	rowNames: A list of site (row) names. If None, inherits from the index of Y.
	scaling: Which scaling should be used. See online documentation.

	Attributes
	---------
	r_w: Row weights
	c_w: Column weights
	evals: Constrained eigenvalues
	U: Constrained eigenvectors
	resid: A pandas.DataFrame of residuals from constrained ordination
	spScores: A pandas.DataFrame of species scores
	siteScores: A pandas.DataFrame of site scores
	siteFitted: A pandas.DataFrame of constrained site scores
	varScores: A pandas.DataFrame containing predictor loadings on each axis
	res_evals: Eigenvalues for unconstrained axes
	res_evecs: Eigenvectors for unconstrained axes

	Methods
	--------
	summary(): 
		Returns a dataframe of summary statistics for each axis
	anova(nperm=999):
		Permutation test of global CCA significance
	triplot(xax=1, yax=2):
		Creates a triplot of components/axes. 
			xax: which component to put on the x-axis
			yax: which component to put on the y-axis

	Example
	--------
	import ecopy as ep

	varespec = ep.load_data('varespec')
	varechem = ep.load_data('varechem')

	cca_fit = ep.cca(varespec, varechem)
	print(cca_fit.summary())
	print(cca_fit.anova())
	cca_fit.triplot()
	"""
	def __init__(self, Y, X, varNames_y=None, varNames_x=None, rowNames=None, scaling=1):
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
		if scaling not in [1,2]:
			msg = 'Scaling must be 1 or 2'
			raise ValueError(msg)
		tolerance = 1E-6
		y_mat = np.array(Y, dtype='float')
		if isinstance(X, np.ndarray):
			x_mat = X
		elif isinstance(X, DataFrame):
			x_mat = np.array(get_dummies(X), dtype='float')
		nrow = X.shape[0]
		ncol_x = X.shape[1]
		ncol_y = Y.shape[1]
		tot = y_mat.sum()
		self.r_w = y_mat.sum(axis=1).reshape(1, nrow) / tot
		self.c_w = y_mat.sum(axis=0).reshape(ncol_y, 1) / tot
		x_mu = self.r_w.dot(x_mat)
		x_mu =  x_mu.reshape(1, ncol_x)
		ones = np.ones(nrow).reshape(nrow, 1)
		mu_mat = ones.dot(x_mu)
		D = x_mat - mu_mat
		sd = np.sqrt(self.r_w.dot(D**2))
		scale_mat = np.diag(1/sd.flatten())
		x_scale = D.dot(scale_mat)
		O_mat = y_mat  / tot
		F_mat = self.r_w.T.dot(self.c_w.T)
		Q_mat = (O_mat - F_mat) / np.sqrt(F_mat)
		W = np.diag(self.r_w.flatten())
		B = np.linalg.pinv(x_scale.T.dot(W).dot(x_scale)).dot(x_scale.T).dot(W**0.5).dot(Q_mat)
		Yhat = (W**0.5).dot(x_scale).dot(B)
		Syy = Yhat.T.dot(Yhat)
		evals, evecs = np.linalg.eig(Syy)
		idx = evals.argsort()[::-1]
		self.evals = np.real(evals[idx])
		self.U = np.real(evecs[:,idx])
		self.evals = self.evals[self.evals>tolerance]
		self.U = self.U[:,self.evals>tolerance]
		Uhat = Q_mat.dot(self.U).dot(np.diag(self.evals**-0.5))
		Wcol = np.diag(self.c_w.flatten()**-0.5)
		Wrow = np.diag(self.r_w.flatten()**-0.5)
		V = Wcol.dot(self.U)
		Vhat = Wrow.dot(Uhat)
		F = Vhat.dot(np.diag(self.evals**0.5))
		Fhat = V.dot(np.diag(self.evals**0.5))
		CAlist = ['CA Axis {0}'.format(i) for i in range(1, len(self.evals)+1)]
		self.resid = np.array(DataFrame(Q_mat - Yhat, columns=varNames_y, index=rowNames))
		self.res_evals, self.res_evecs = np.linalg.eig(self.resid.T.dot(self.resid))
		self.res_evals = np.real(self.res_evals[self.res_evals>1E-9])
		self.y_mat = y_mat
		self.x_mat = x_mat
		if scaling==1:
			Wrow = np.diag(self.r_w.flatten()**0.5)
			self.spScores = DataFrame(V, columns=CAlist, index=varNames_y)
			self.siteScores = DataFrame(F, columns=CAlist, index=rowNames)
			self.siteFitted = DataFrame(Wrow.dot(Yhat).dot(self.U), columns=CAlist, index=rowNames)
			mu_z =  self.r_w.dot(self.siteFitted)
			mu_z_mat = ones.dot(mu_z)
			D_z = self.siteFitted - mu_z_mat
			sd_z = np.sqrt(self.r_w.dot(D_z**2))
			scale_z_mat = np.diag(1/sd_z.flatten())
			scaled_Z= D_z.dot(scale_z_mat)
			self.varScores = DataFrame(x_scale.T.dot(Wrow).dot(scaled_Z).dot(np.diag(self.evals**0.5)), columns=CAlist, index=varNames_x)
		if scaling==2:
			Wrow = np.diag(self.r_w.flatten()**0.5)
			self.spScores = DataFrame(Fhat, columns=CAlist, index=varNames_y)
			self.siteScores = DataFrame(Vhat, columns=CAlist, index=rowNames)
			self.siteFitted = DataFrame(Wrow.dot(Yhat).dot(self.U).dot(np.diag(self.evals**-0.5)), columns=CAlist, index=rowNames)
			mu_z =  self.r_w.dot(self.siteFitted)
			mu_z_mat = ones.dot(mu_z)
			D_z = self.siteFitted - mu_z_mat
			sd_z = np.sqrt(self.r_w.dot(D_z**2))
			scale_z_mat = np.diag(1/sd_z.flatten())
			scaled_Z= D_z.dot(scale_z_mat)
			self.varScores = DataFrame(x_scale.T.dot(Wrow).dot(scaled_Z), columns=CAlist, index=varNames_x)

	def summary(self):
		print('Constrained variance = {0:.3}'.format(np.sum(self.evals)))
		print('Unconstrained varience = {0:.3}'.format(np.sum(self.res_evals)))
		names = ['CCA {0}'.format(x) for x in range(1, len(self.evals)+1)]
		data = np.vstack((np.round(self.evals, 3), np.round(self.evals/self.evals.sum(),3)))
		SumTable1 = DataFrame(data, index = ['Variance', 'Prop. Variance'], columns=names)
		print('Constrained Axes')
		print(SumTable1)
		print('\n')
		names2 = ['CA {0}'.format(x) for x in range(1, len(self.res_evals)+1)]
		data2 = np.vstack((np.round(self.res_evals, 3), np.round(self.res_evals/self.res_evals.sum(),3)))
		SumTable2 = DataFrame(data2, index = ['Variance', 'Prop. Variance'], columns=names2)
		print('Unconstrained Axes')
		print(SumTable2)

	def anova(self, nperm=999):
		constrained = np.sum(self.evals)
		unconstrained = np.sum(self.res_evals)
		Fobs = (constrained/len(self.evals)) / (unconstrained/len(self.res_evals))
		Fperm = np.empty(nperm)
		for i in range(nperm):
			n = self.y_mat.shape[0]
			nrow = n
			ncol_y = self.y_mat.shape[1]
			ncol_x = self.x_mat.shape[1]
			idx = np.random.choice(n, n, replace=False)
			y_perm = self.y_mat[idx,:]
			tot = y_perm.sum()
			r_w = y_perm.sum(axis=1).reshape(1, nrow) / tot
			c_w = y_perm.sum(axis=0).reshape(ncol_y, 1) / tot
			x_mu = r_w.dot(self.x_mat)
			x_mu =  x_mu.reshape(1, ncol_x)
			ones = np.ones(nrow).reshape(nrow, 1)
			mu_mat = ones.dot(x_mu)
			D = self.x_mat - mu_mat
			sd = np.sqrt(r_w.dot(D**2))
			scale_mat = np.diag(1/sd.flatten())
			x_scale = D.dot(scale_mat)
			O_mat = y_perm  / tot
			F_mat = r_w.T.dot(c_w.T)
			Q_mat = (O_mat - F_mat) / np.sqrt(F_mat)
			W = np.diag(r_w.flatten())
			B = np.linalg.pinv(x_scale.T.dot(W).dot(x_scale)).dot(x_scale.T).dot(W**0.5).dot(Q_mat)
			Yhat = (W**0.5).dot(x_scale).dot(B)
			Syy = Yhat.T.dot(Yhat)
			evals_perm, evecs = np.linalg.eig(Syy)
			evals_perm = np.real(evals_perm[evals_perm.argsort()[::-1]])
			evals_perm = evals_perm[evals_perm>1E-6]
			resid = np.array(DataFrame(Q_mat - Yhat))
			res_evals, res_evecs = np.linalg.eig(resid.T.dot(resid))
			res_evals = np.real(res_evals[res_evals>1E-6])
			c_perm = np.sum(evals_perm)
			u_perm = np.sum(res_evals)
			Fperm[i] = (c_perm/len(evals_perm)) / (u_perm/len(res_evals))
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
		for i in range(self.varScores.shape[0]):
			ax.arrow(0, 0, self.varScores.iloc[i, xplot], self.varScores.iloc[i, yplot], color='b', head_width=0.1)
			ax.text(self.varScores.iloc[i,xplot]*1.2, self.varScores.iloc[i,yplot]*1.2, self.varScores.index.values[i], color='b', ha='center', va='center')
		xmins = (self.spScores.iloc[:,xplot].min(), self.siteScores.iloc[:,xplot].min(), self.varScores.iloc[:,xplot].min()*1.2)
		xmax = (self.spScores.iloc[:,xplot].max(), self.siteScores.iloc[:,xplot].max(), self.varScores.iloc[:,xplot].max()*1.2)
		ymins = (self.spScores.iloc[:,yplot].min(), self.siteScores.iloc[:,yplot].min(), self.varScores.iloc[:,yplot].min()*1.2)
		ymax = (self.spScores.iloc[:,yplot].max(), self.siteScores.iloc[:,yplot].max(), self.varScores.iloc[:,yplot].max()*1.2)
		ax.set_xlabel('CA {0}'.format(xax))
		ax.set_ylabel('CA {0}'.format(yax))
		ax.set_xlim([min(xmins), max(xmax)])
		ax.set_ylim([min(ymins), max(ymax)])
		plt.show()
		




