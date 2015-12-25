import numpy as np
from pandas import DataFrame, get_dummies
import matplotlib.pyplot as plt
from ..base_funcs import wt_scale

class hillsmith(object):
	'''
	Docstring for function ecopy.hillsmith
	====================
	Conducts ordination on a matrix of mixed data types (both quantitative and qualitative)
		following Hill and Smith (1976). If DataFrame is all quantitative, this is equivalent to
		principle components analysis. If DataFrame is all qualitative, this is equivalent 
		to multiple correspondance analysis.

	Use
	----
	hillsmith(mat, wt_r=None, ndim=2)

	Returns an object of class hillsmith

	Parameters
	----------
	mat:  Data for ordination. Must be a pandas DataFrame. 
		Observations as rows and descriptors as columns.
		Rows must be > columns.
	wt_r: Optional vector of row weights. Can be a numpy.ndarray or 
		a list.
	ndim: Number of axes and components to extract

	Attributes
	---------
	pr_axes: DataFrame of principle axis loadings (columns)
	row_coords: DataFrame of row coordinates on each principle axis
	pr_components: DataFrame of principle component loadings (rows)
	column_coords: DataFrame of column coordinates on each component
	evals: eigenvalues for each axis

	Methods
	--------
	summary(): Returns a dataframe of summary 
		statistics for each axis
	biplot(invert=False, xax=1, yax=2, obsNames=True):
		creates a biplot of components/axes 
			invert: If False, the plots rows as points and columns
				as arrows. If True, plots columns as points
				and rows as arrows.
			xax: which component to put on the x-axis
			yax: which component to put on the y-axis
			obsNames: tells whether to plot the names of each observation (True by default)

	Example
	--------
	import ecopy as ep
	dune_env = ep.load_data('dune_env')
	dune_env = dune_env[['A1', 'Moisture', 'Manure', 'Use', 'Management']]
	print(ep.hillsmith(dune_env).summary().iloc[:,:2])
	ep.hillsmith(dune_env).biplot(obsNames=False, invert=False)
	'''
	def __init__(self, mat, wt_r=None, ndim=2):

		def corfunc(x):
			column = x[0]
			component = x[1]
			if columnType[column]=='q':
				w = components[:,component] * wt_r * modMat[:,columnIndex==column].ravel()
				return np.sum(w)**2
			else:
				x = components[:,component]*wt_r
				qual = np.array(df.iloc[:,column])
				groups = np.unique(qual)
				denom = []
				num = []
				for g in groups:
					denom.append(wt_r[qual==g].sum())
					num.append(x[qual==g].sum())
				denom = np.array(denom)
				num = np.array(num)
				div = num / denom
				return np.sum(denom*div*div)

		def ordfunc(mat, wt_col, wt_row, ndim=2):
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

		if not isinstance(mat, DataFrame):
			msg = 'Data must be pandas.DataFrame'
			raise ValueError(msg)
		if isinstance(mat, DataFrame):
			if mat.isnull().any().any():
				msg = 'DataFrame contains null values'
				raise ValueError(msg)
		df = mat
		nr = df.shape[0]
		nc = df.shape[1]
		if nr < nc:
			msg = 'Number of columns cannot exceed number of columns'
			raise ValueError(msg)
		columnType = np.array([""]*nc)
		for column in range(nc):
			w1 = 'q'
			if df.iloc[:,column].dtype=='int':
				df.iloc[:,column] = df.iloc[:,column].apply(float)
			elif df.iloc[:,column].dtype=='object':
				w1 = 'f'
			columnType[column] = w1
		modMat = np.array([0]*nr).reshape(nr, 1)
		modNames = ['0']
		if wt_r is None:
			wt_r = np.array([1.]*nr)
		wt_r = np.array(wt_r) / np.array(wt_r).sum()
		wt_c = np.array(np.nan)
		columnIndex = np.array(np.nan)
		columnID = 0
		for column in range(nc):
			if columnType[column]=='q':
				modMat = np.hstack((modMat, wt_scale(df.iloc[:,column], wt=wt_r, bias=1).reshape(nr, 1)))
				modNames.append(df.columns[column])
				wt_c = np.append(wt_c, 1)
				columnIndex = np.append(columnIndex, columnID)
				columnID += 1
			elif columnType[column]=='f':
				colName = df.columns[column]
				w = get_dummies(df.iloc[:,column])
				levels = [colName + ':' + x for x in w.columns]
				temp_wt = wt_r.dot(w)
				w = w.apply(lambda x: x/temp_wt - 1, axis=1)
				wt_c = np.append(wt_c, temp_wt)
				modMat = np.hstack((modMat, w))
				modNames.extend(levels)
				columnIndex = np.append(columnIndex, np.array([columnID]*len(levels)))
				columnID += 1
		modMat = modMat[:,1:]
		wt_c = wt_c[1:]
		columnIndex = columnIndex[1:]
		modNames = modNames[1:]
		idMat = np.zeros((nc, ndim, 2), dtype='int')
		idMat[:,:,0] = np.repeat(np.arange(nc).reshape(nc, 1), ndim, axis=1)
		idMat[:,:,1] = np.repeat(np.arange(ndim).reshape(1, ndim), nc, axis=0)
		axes, rowcoords, colcoords, components, evals = ordfunc(modMat, wt_c, wt_r)
		corrMat = DataFrame(np.apply_along_axis(corfunc, 2, idMat), index=df.columns)
		corrMat.columns = ['Comp1', 'Comp2']
		self.pr_axes = DataFrame(axes, index=modNames)
		axis_names = ['Axis {0}'.format(str(i+1)) for i in range(axes.shape[1])]
		self.pr_axes.columns = axis_names
		self.row_coords = DataFrame(rowcoords, index=df.index)
		self.row_coords.columns = axis_names
		self.pr_components = DataFrame(components, index=df.index)
		component_names = ['Component {0}'.format(str(i+1)) for i in range(components.shape[1])]
		self.pr_components.columns = component_names
		self.column_coords = DataFrame(colcoords, index=modNames)
		self.column_coords.columns = component_names
		self.evals = evals

	def summary(self):
		axis_names = ['Axis {0}'.format(str(i+1)) for i in range(len(self.evals))]
		sds = np.sqrt(self.evals)
		props = self.evals / self.evals.sum()
		cumprop = np.cumsum(self.evals) / self.evals.sum()
		summDF = DataFrame(np.vstack((sds, props, cumprop)), index=['Std. Dev', 'Prop Var', 'Cum Var'])
		summDF.columns = axis_names
		return summDF

	def biplot(self, invert=False, xax=1, yax=2, obsNames=True):
		points = self.row_coords
		arrows = self.pr_axes
		if invert:
			points = self.column_coords
			arrows = self.pr_components
		f, ax = plt.subplots()
		ax.axvline(0, ls='solid', c='k')
		ax.axhline(0, ls='solid', c='k')
		if obsNames:
			ax.scatter(points.iloc[:,xax-1], points.iloc[:,yax-1], s=0)
			for i in range(points.shape[0]):
				plt.text(points.iloc[i,xax-1], points.iloc[i,yax-1], points.index.values[i], ha = 'center', va = 'center')
		else:
			ax.scatter(points.iloc[:,xax-1], points.iloc[:,yax-1])
		for i in range(arrows.shape[0]):
			ax.arrow(0, 0, arrows.iloc[i,xax-1], arrows.iloc[i,yax-1], color = 'red', head_width=.05)
			ax.text(arrows.iloc[i, xax-1]*1.2, arrows.iloc[i,yax-1]*1.2, arrows.index.values[i], color = 'red', ha = 'center', va = 'center')
		xmax = max(np.amax(points.iloc[:,xax-1]), np.amax(arrows.iloc[:,xax-1]))
		xmin = min(np.min(points.iloc[:,xax-1]), np.min(arrows.iloc[:,xax-1]))
		ymax = max(np.amax(points.iloc[:,yax-1]), np.amax(arrows.iloc[:,yax-1]))
		ymin = min(np.amin(points.iloc[:,yax-1]), np.amin(arrows.iloc[:,yax-1]))
		ax.set_xlim([xmin + 0.15*xmin, xmax+0.15*xmax])
		ax.set_ylim([ymin + 0.15*ymin, ymax+0.15*ymax])
		ax.set_xlabel('Axis {!s}'.format(xax))
		ax.set_ylabel('Axis {!s}'.format(yax))
		plt.show()


