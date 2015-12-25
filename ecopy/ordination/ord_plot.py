from scipy.spatial import ConvexHull
import numpy as np
from pandas import DataFrame, Series
import matplotlib.pyplot as plt

def ord_plot(x, groups, y=None, colors=None, type='Hull', label=True, showPoints=True, xlab='Axis 1', ylab='Axis 2'):
	'''
	Docstring for function ecopy.ord_plot
	====================
	Delineates different groups in ordination (or regular) space.

	Use
	----
	ord_plot(x, groups, y=None, colors=None, type='Hull', label=True, showPoints=True, xlab='Axis 1', ylab='Axis 2')

	Parameters
	----------
	x: a numpy.ndarray, pandas.DataFrame, or pandas.Series containing coordinates to be plotted. Can be either a one
		or two column matrix. If only one column, then y must be specified.
	groups: a list, pandas.DataFrame, or pandas.Series containing a factor denoting group identification
	y: a numpy.ndarray, pandas.DataFrame, or pandas.Series containing the y coordinates to be plotted. Can
		only have one column, and must be specified is x is only one column.
	colors: a string, list, pandas.Series, or pandas.DataFrame giving custom colors for each group. Otherwise
		default colors are used.
	type: either 'Hull' or 'Line'. 'Hull' produces a convex hull, whereas 'Line' produces lines connected to the centroid for each
		point.
	label: Whether or not a label should be shown at the center of each group.
	showPoints: Whether or not the points should be shown.
	xlab: Label for the x-axis.
	ylab: Label for the y-axis.

	Example
	--------
	import numpy as np
	import ecopy as ep
	import matplotlib.pyplot as plt

	nObs = 10
	X = np.random.normal(0, 1, 10*2)
	Y = np.random.normal(0, 1, 10*2)
	GroupID = ['A']*nObs + ['B']*nObs

	Z = np.vstack((X, Y)).T

	# x contains two columns. convex hull. custom colors.
	ep.ord_plot(x=Z, groups=GroupID, colors=['r', 'b'])

	# x and y coordinates in different matrices. line plot.
	ep.ord_plot(x=X, y=Y, groups=GroupID, type='Line', xlab='PC1', ylab='PC2', showPoints=False, label=False)
	'''
	if not isinstance(x, (np.ndarray, DataFrame, Series)):
		msg = 'x must be a numpy.ndarray, pandas.DataFrame, or pandas.Series'
		raise ValueError (msg)
	x = np.array(x)
	if y is not None:
		if not isinstance(y, (np.ndarray, DataFrame, Series)):
			msg = 'y must be a numpy.ndarray, pandas.DataFrame, or pandas.Series'
		z = np.vstack((x, y)).T
		if z.shape[1]!=2:
			msg = 'x and y can only have one column if y is provided'
			raise ValueError(msg)
	else:
		if x.shape[1] != 2:
			msg = 'x must be a n x 2 matrix'
			raise ValueError(msg)
		z = x
	if type not in ['Hull', 'Line']:
		msg = 'type must either be Hull or Line'
		raise ValueError(msg)
	if not isinstance(groups, (list, DataFrame, Series)):
		msg = 'groups must be either a list, pandas.DataFrame, or pandas.Series'
		raise ValueError(msg)
	g = list(groups)
	unique_g = list(set(g))
	ng = len(unique_g)
	if colors is None:
		cmap = plt.get_cmap('Paired')
		cID = [cmap(i) for i in np.linspace(0, 1, ng)]
	else:		
		if not isinstance(colors, (list, Series, DataFrame)):
			msg = 'colors must be a string, list, pandas.DataFrame, or pandas.Series'
		if isinstance(colors, (list, Series, DataFrame)):
			cID = list(colors)
		if isinstance(colors, str):
			cID = [colors]*ng
	if len(cID) != len(unique_g):
		msg = 'list of colors must equal number of groups'
		raise ValueError(msg)
	if type=='Hull':	
		f, ax = plt.subplots()
		for j in range(len(unique_g)):
			tempBOOL = np.array(groups)==unique_g[j]
			tempX = z[tempBOOL,:]
			tempHull = ConvexHull(tempX)
			for simplex in tempHull.simplices:
				ax.plot(tempX[simplex, 0], tempX[simplex, 1], '-', c=cID[j], lw=1)
			if showPoints:
				ax.plot(tempX[:,0], tempX[:,1], 'o', mfc=cID[j], mew=1, label=unique_g[j])
			if label:
				Xcent = tempX[:,0].mean()
				Ycent = tempX[:,1].mean()
				ax.text(Xcent, Ycent, unique_g[j], color=cID[j], ha='center', va='center', bbox={'facecolor': 'white', 'pad': 10, 'edgecolor': cID[j], 'lw': 1})
		ax.set_xlabel(xlab)
		ax.set_ylabel(ylab)	
		plt.show()
	if type=='Line':
		f, ax = plt.subplots()
		for j in range(len(unique_g)):
			tempBOOL = np.array(groups)==unique_g[j]
			tempX = z[tempBOOL,:]
			Xcent = tempX[:,0].mean()
			Ycent = tempX[:,1].mean()
			for l in range(tempX.shape[0]):
				ax.plot([Xcent, tempX[l,0]], [Ycent, tempX[l,1]], c=cID[j], lw=1)
			if showPoints:
				ax.plot(tempX[:,0], tempX[:,1], 'o', mfc=cID[j], mew=1, label=unique_g[j])
			if label:
				ax.text(Xcent, Ycent, unique_g[j], color=cID[j], ha='center', va='center', bbox={'facecolor': 'white', 'pad': 10, 'edgecolor': cID[j], 'lw': 1})
		ax.set_xlabel(xlab)
		ax.set_ylabel(ylab)
		plt.show()
