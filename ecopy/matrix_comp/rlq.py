import numpy as np
from pandas import DataFrame, get_dummies
import matplotlib.pyplot as plt
from ..base_funcs import wt_scale

class rlq(object):
	"""
	Docstring for function ecopy.rlq
	====================
	Conducts RLQ analysis for an environmental matrix (R),
		a species matrix (L), and a trait matrix (Q)

	Use
	----
	rlq(R, L, Q, ndim=2)

	Returns an object of class rlq

	Parameters
	----------
	R: A pandas.DataFrame containing environmental data (site x environmental data)
	L:  A numpy.ndarray or pandas.DataFrame containing species presence/absence
		or abundances at each site (site x species).
	Q: A pandas.DaraFrame containing species traits (species x trait)
	ndim: Number of axes and components to save

	Attributes
	---------
	traitVecs: DataFrame of trait loadings
	normedTraits: DataFrame of species (traits) scores, normalized to unit length
	envVecs: DataFrame of environmental loadings
	normedEnv: DataFrame of site (environment) scores, normalized to unit length
	evals: eigenvalues for each axis

	Methods
	--------
	summary(): Returns a dataframe of summary 
		statistics for each axis
	biplot(xax=1, yax=2):
		creates a biplot of components/axes. The output
		is a four-panel graph showing species loadings,
		site loadings, trait vectors, and environmental 
		vectors on different panels
			xax: which component to put on the x-axis
			yax: which component to put on the y-axis

	Example
	--------
	import ecopy as ep

	traits = ep.load_data('avi_traits')
	env = ep.load_data('avi_env')
	sp = ep.load_data('avi_sp')

	rlq_output = ep.rlq(env, sp, traits)
	print(rlq_ouput.summary().iloc[:,:3])
	rlq_output.biplot()
	"""
	def __init__(self, R, L, Q, ndim=2):
		if not isinstance(R,DataFrame):
			msg = 'Matrix R must be a pandas.DataFrame'
			raise ValueError(msg)
		if not isinstance(L, DataFrame):
			msg = 'Matrix L must be a pandas.DataFrame'
			raise ValueError(msg)
		if L.dtypes.any()=='object':
			msg ='Matrix L must be only numeric'
			raise ValueError(msg)
		if not isinstance(Q, DataFrame):
			msg = 'Matrix Q must be a pandas.DataFrame'

		nr_sp = L.shape[0]
		nc_sp = L.shape[1]
		nr_e = R.shape[0]
		nc_e = R.shape[1]
		nr_t = Q.shape[0]
		nc_t = Q.shape[1]

		if nr_e != nr_sp:
			msg = 'Matrices R and L must have the same number of rows'
			raise ValueError(msg)
		if nr_t != nc_sp:
			msg = 'The number of matrix L columns must equal  the number of matrix Q rows'
		if nr_e < nc_e:
			msg = 'Number of columns cannot exceed number of rows in environment matrix'
			raise ValueError(msg)
		if nr_t < nc_t:
			msg = 'Number of columns cannot exceed number of rows in trait matrix'
			raise ValueError(msg)

		Lmat = np.array(L, dtype='float')
		if Lmat.sum(axis=0).any()==0:
			msg = 'Matrix L has at least one empty column'
			raise ValueError(msg)
		if Lmat.sum(axis=1).any()==0:
			msg = 'Matrix L has at least one empty row'
			raise ValueError(msg)

		Lmat = Lmat / Lmat.sum()
		row_w = Lmat.sum(axis=1)
		col_w = Lmat.sum(axis=0)
		Lmat = np.apply_along_axis(lambda x: x / row_w, 0, Lmat)
		Lmat = np.apply_along_axis(lambda x: x/col_w, 1, Lmat) - 1
		envMat, envNames, envWeights = dummyMat(R, nr_e, nc_e, row_w)
		traitMat, traitNames, traitWeights = dummyMat(Q, nr_t, nc_t, col_w)
		RLQmat = envMat.T.dot(np.diag(row_w)).dot(Lmat).dot(np.diag(col_w)).dot(traitMat)
		RLQmat = DataFrame(RLQmat, index=envNames)
		RLQmat.columns = traitNames
		axes, rowcoords, colcoords, components, self.evals = ordfunc(RLQmat, traitWeights, envWeights, ndim)
		self.evals = np.real(self.evals)
		self.traitVecs = np.apply_along_axis(lambda x: x*traitWeights, 0, axes)
		traitScores = DataFrame(np.real(traitMat.dot(self.traitVecs)), index=L.columns)
		traitScores.columns = ['Trait Axis {0}'.format(x) for x in range(1,ndim+1)]
		self.traitVecs = DataFrame(np.real(axes), index=traitNames)
		self.traitVecs.columns = ['Trait Vector {0}'.format(x) for x in range(1,ndim+1)]
		self.normedTraits = normalize(traitScores, col_w)
		self.envVecs = np.apply_along_axis(lambda x: x*envWeights, 0, components)
		envScores = DataFrame(np.real(envMat.dot(self.envVecs)), index=L.index)
		envScores.columns = ['Environment Component {0}'.format(x) for x in range(1,ndim+1)]
		self.normedEnv = normalize(envScores, row_w)
		self.envVecs = DataFrame(np.real(components), index=envNames)
		self.envVecs.columns = ['Environmental Vector {0}'.format(x) for x in range(1,ndim+1)]

	def summary(self):
		axis_names = ['Axis {0}'.format(str(i+1)) for i in range(len(self.evals))]
		sds = np.sqrt(self.evals)
		props = self.evals / self.evals.sum()
		cumprop = np.cumsum(self.evals) / self.evals.sum()
		summDF = DataFrame(np.vstack((sds, props, cumprop)), index=['Std. Dev', 'Prop Var', 'Cum Var'])
		summDF.columns = axis_names
		return summDF

	def biplot(self, xax=1, yax=2):
		f, ax = plt.subplots(2, 2, figsize=(12, 12))
		ax[0,0].plot(self.normedTraits.iloc[:,xax-1], self.normedTraits.iloc[:,yax-1],'o', ms=0)
		for i in range(self.normedTraits.shape[0]):
			ax[0,0].text(self.normedTraits.iloc[i, xax-1], self.normedTraits.iloc[i, yax-1], self.normedTraits.index.values[i], ha='center', va='center')
			ax[0,0].set_xlabel(self.normedTraits.columns[0])
			ax[0,0].set_ylabel(self.normedTraits.columns[1])
		ax[0,1].plot(self.normedEnv.iloc[:,xax-1], self.normedEnv.iloc[:,yax-1],'o', ms=0)
		for i in range(self.normedEnv.shape[0]):
			ax[0,1].text(self.normedEnv.iloc[i, xax-1], self.normedEnv.iloc[i, yax-1], self.normedEnv.index.values[i], ha='center', va='center')
			ax[0,1].set_xlabel(self.normedEnv.columns[0])
			ax[0,1].set_ylabel(self.normedEnv.columns[1])
		for i in range(self.traitVecs.shape[0]):
			xlim = [np.min(self.traitVecs.iloc[:,xax-1])*1.1, np.max(self.traitVecs.iloc[:,xax-1])*1.1]
			ylim = [np.min(self.traitVecs.iloc[:,yax-1])*1.1, np.max(self.traitVecs.iloc[:,yax-1])*1.1]
			ax[1,0].arrow(0, 0, self.traitVecs.iloc[i,xax-1], self.traitVecs.iloc[i,yax-1], color='red', head_width=np.ptp(xlim)*0.01)
			ax[1,0].text(self.traitVecs.iloc[i,xax-1]*1.1, self.traitVecs.iloc[i,yax-1]*1.1, self.traitVecs.index.values[i], color='r', ha='center', va='center')
			ax[1,0].set_xlim(xlim)
			ax[1,0].set_ylim(ylim)
			ax[1,0].set_xlabel(self.traitVecs.columns[0])
			ax[1,0].set_ylabel(self.traitVecs.columns[1])
		for i in range(self.envVecs.shape[0]):
			xlim = [np.min(self.envVecs.iloc[:,xax-1])*1.1, np.max(self.envVecs.iloc[:,xax-1])*1.1]
			ylim = [np.min(self.envVecs.iloc[:,yax-1])*1.1, np.max(self.envVecs.iloc[:,yax-1])*1.1]
			ax[1,1].arrow(0, 0, self.envVecs.iloc[i,xax-1], self.envVecs.iloc[i,yax-1], color='blue', head_width=np.ptp(xlim)*0.01)
			ax[1,1].text(self.envVecs.iloc[i,xax-1]*1.1, self.envVecs.iloc[i,yax-1]*1.1, self.envVecs.index.values[i], color='blue', ha='center', va='center')
			ax[1,1].set_xlim(xlim)
			ax[1,1].set_ylim(ylim)
			ax[1,1].set_xlabel(self.envVecs.columns[0])
			ax[1,1].set_ylabel(self.envVecs.columns[1])
		plt.subplots_adjust(wspace=0.3, hspace=0.3)
		plt.show()

def dummyMat(matrix, rows, columns, row_w):
	columnType = np.array([""]*columns)
	weights = np.array(np.nan)
	for col in range(columns):
		id1 = 'q'
		if matrix.iloc[:,col].dtype=='int':
			matrix.iloc[:,col] = matrix.iloc[:,col].apply(float)
		elif matrix.iloc[:,col].dtype=='object':
			id1 = 'f'
		columnType[col] = id1
	modMat = np.array([0]*rows).reshape(rows, 1)
	modNames = ['0']
	for col in range(columns):
		if columnType[col]=='q':
			modMat = np.hstack((modMat, wt_scale(matrix.iloc[:,col], row_w, bias=1).reshape(rows, 1)))
			modNames.append(matrix.columns[col])
			weights = np.append(weights, 1)
		elif columnType[col]=='f':
			colName = matrix.columns[col]
			w = get_dummies(matrix.iloc[:,col])
			levels = [colName + ':' + x for x in w.columns]
			temp_wt = row_w.dot(w)
			w = w.apply(lambda x: x/temp_wt - 1, axis=1)
			modMat = np.hstack((modMat, w))
			modNames.extend(levels)
			weights = np.append(weights, temp_wt)
	return modMat[:,1:], modNames[1:], weights[1:]

def ordfunc(mat, wt_col, wt_row, ndim):
			X = np.array(mat)
			X2 = np.apply_along_axis(lambda x: x*np.sqrt(wt_row), 0, X)
			X2 = np.apply_along_axis(lambda x: x*np.sqrt(wt_col), 1, X2)
			X2 = X2.T.dot(X2)
			evals, evecs = np.linalg.eig(X2)
			evals = evals[evals.argsort()[::-1]]
			evecs = evecs[:,evals.argsort()[::-1]]
			evecs = evecs[:,:ndim]
			sds = np.sqrt(evals)[:ndim]
			axiswt = 1/np.sqrt(wt_col)
			axes = np.apply_along_axis(lambda x: x*axiswt, 0, evecs)
			rowcoords = np.apply_along_axis(lambda x: x*wt_col, 1, X)
			rowcoords = rowcoords.dot(axes)
			colcoords = np.apply_along_axis(lambda x: x*sds, 1, axes)
			components = np.apply_along_axis(lambda x: x/sds, 1, rowcoords)
			return axes, rowcoords, colcoords, components, evals

def normalize(X, w):
	Z = np.array(X)
	norms = np.apply_along_axis(lambda x: np.sqrt(np.sum(x*x*w) / np.sum(w)), 0, Z)
	normedMat = np.apply_along_axis(lambda x: x / norms, 1, Z)
	normedMat = DataFrame(normedMat, index=X.index)
	normedMat.columns = X.columns
	return normedMat