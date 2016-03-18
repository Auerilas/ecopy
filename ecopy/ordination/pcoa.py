import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as py
from warnings import warn

class pcoa(object):
	'''
	Docstring for function ecopy.pcoa
	====================
	Conducts principle coordinate analysis of a user-supplied
		distance matrix. Matrix must be square symmetric

	Use
	----
	pcoa(x, correction=None, siteNames=None)

	Returns an object of class pcoa

	Parameters
	----------
	x:  Data for ordination. Should be a pandas.DataFrame or
		numpy.ndarray square, symmetric matrix of distances.
		See ecopy.distance for calculation of distances.
		See online documentation for details about the
		method
	correction: Which correction for negative eigenvalues
		should be applied. Accepts either '1' or '2'. 
		See online documentation for details.
	siteNames: A list of site names used in biplots

	Attributes
	---------
	evals: Eigenvalues of the transformed distance matrix
	U: Eigenvectors of the transformed distance matrix.
		Each eigenvector has already been multiplied
		by the square root of its eigenvalue
	correction: The correction applied to correct for negative
		eigenvalues
	
	Methods
	--------
	summary(): Returns a pandas.DataFrame summary
		table of PCoA axes
	biplot(coords=False, xax=1, yax=2, descriptors=None, descripNames=None):
	   	Produces a biplot of given PCoA axes
	   	coords: If true, return coordinates for objects and descriptors for
	   		custom plotting
		xax: An integer specifying which PCoA axis should appear
			on the x-axis (Defaults to 1)
		yax: An integer specifying which PCoA axis should appear
			on the y-axis (Defaults to 2)
		descriptors: A pandas.DataFrame or numpy.ndarray
			of descriptors to project onto the biplot. Descriptors
			are all standardized prior to projection
		descripNames: A list of names for the descriptors. If
			not specified, will be taken from the column names
			if descriptors is a pandas.DataFrame
	shepard(xax=1, yax=2): Produces a shepard diagram
		plotting the euclidean distances in reduced space
		against observed distances from the original matrix
	
	Example
	--------
	import ecopy as ep

	BCI = ep.load_data('BCI')
	brayD = ep.distance(BCI, method='bray', transform='sqrt')
	pc1 = ep.pcoa(brayD)
	print(pc1.summary())
	pc1.biplot()
	pc1.shepard()
	'''
	def __init__(self, x, correction=None, siteNames=None):
		if not isinstance(x, (DataFrame, np.ndarray)):
			msg = 'Data must either be pandas.DataFrame or nump.ndarray'
			raise ValueError(msg)
		if isinstance(x, DataFrame):
			if (x.dtypes == 'object').any():
				msg = 'DataFrame can only contain numeric values'
				raise ValueError(msg)
			y = np.array(x)
		if isinstance(x, np.ndarray):
			y = x
		if np.isnan(y).any():
			msg = 'Distance matrix contains null values'
			raise ValueError(msg)
		if y.any() < 0:
			msg ='Distance matrix cannot contain negative values'
			raise ValueError(msg)
		if y.shape[0] != y.shape[1]:
			msg = 'Distance matrix must be square'
			raise ValueError(msg)
		if not np.allclose(y.T, y):
			msg ='Distance matrix must be symmetric'
			raise ValueError(msg)
		A = -0.5*np.square(y.astype('float'))
		n = float(y.shape[0])
		ones = np.ones(n)[np.newaxis].T
		I = np.eye(n)
		D = (I - ones.dot(ones.T)/n).dot(A).dot(I - ones.dot(ones.T)/n)
		self.evals, self.U = np.linalg.eig(D)
		self.evals = self.evals.real
		self.U = self.U.real
		idx = self.evals.argsort()[::-1]
		self.U = self.U[:,idx]
		if correction is not None:
			if correction not in ['1', '2']:
				msg = "correction must be either '1' or '2'"
				raise ValueError(msg)
		if correction is '1':
			negEvl = np.abs(np.min(self.evals[self.evals < 0]))
			A = -0.5*np.square(y) - negEvl
			np.fill_diagonal(A, 0)
			D = (I - ones.dot(ones.T)/n).dot(A).dot(I - ones.dot(ones.T)/n)
			self.evals, self.U = np.linalg.eig(D)
			idx = self.evals.argsort()[::-1]
			self.U = self.U[:,idx]
			self.correction = negEvl
		if correction is '2':
			mat0 = np.zeros((n,n))
			matI = -1.*np.eye(n)
			d1 = 2.*D
			d2 = -0.5*y
			d2 = -4.*(I - ones.dot(ones.T)/n).dot(d2).dot(I - ones.dot(ones.T)/n)
			t1 = np.concatenate((mat0, matI), axis=0)
			t2 = np.concatenate((d1, d2), axis=0)
			specMat = np.concatenate((t1, t2), axis=1)
			t_evals, t_evecs =np.linalg.eig(specMat)
			posEvl = np.max(np.real(t_evals))
			A = -0.5*(np.square(y.astype('float') + posEvl))
			np.fill_diagonal(A, 0)
			D = (I - ones.dot(ones.T)/n).dot(A).dot(I - ones.dot(ones.T)/n)
			self.evals, self.U = np.linalg.eig(D)
			idx = self.evals.argsort()[::-1]
			self.U = self.U[:,idx]
			self.correction = posEvl
		self.evals = np.round(self.evals[idx], 4)
		self.U = np.round(self.U.dot(np.diag(np.sqrt(self.evals))), 4)
		self.siteLabs = ['Site ' + str(x) for x in range(1, y.shape[0]+1)]
		if isinstance(x, DataFrame):
			self.siteLabs = x.index
		if siteNames is not None:
			self.siteLabs = siteNames
		self.y2 = y

	def summary(self):
		sds = np.sqrt(self.evals)
		props = self.evals / np.sum(self.evals)
		cumSums = np.cumsum(self.evals) / np.sum(self.evals)
		colNames = ['PCoA Axis ' + str(x) for x in range(1, len(self.evals)+1)]
		sumTable = DataFrame(np.vstack((sds, props, cumSums)), index=['Std. Dev', 'Prop.', 'Cum. Prop.'])
		sumTable.columns = colNames
		return sumTable
		
	def biplot(self, coords=False, xax=1, yax=2, descriptors=None, descripNames=None, spCol='r', siteCol='k', spSize=12, siteSize=12):
		if descriptors is not None:
				warn('\nWarning: Descriptors must not be binary.\nIgnore this message if all descriptors are quantitative\n')
				if not isinstance(descriptors, (DataFrame, np.ndarray)):
					msg = 'descriptors must be a pandas.DataFrame or numpy.ndarray'
					raise ValueError(msg)				
				if isinstance(descriptors, DataFrame):
					if (descriptors.dtypes == 'object').any():
						msg = 'DataFrame can only contain numeric values'
						raise ValueError(msg)
				d2 = np.array(descriptors)
				d2 = np.apply_along_axis(lambda x: (x - np.mean(x))/np.std(x, ddof=1), 0, d2)
				if isinstance(descriptors, DataFrame):
					dLabs = descriptors.columns.values
				elif descripNames is not None:
					if len(descripNames) != d2.shape[1]:
						msg = 'descripNames must be equal to the number of columns in descriptors'
						raise ValueError(msg)
					dLabs = descripNames
				else:
					dLabs = ['D'+str(x) for x in range(1, d2.shape[1] + 1)]
				U2 = self.U[:,[xax-1, yax-1]].astype('float')
				U2 = np.apply_along_axis(lambda x: (x - np.mean(x))/np.std(x, ddof=1), 0, U2)
				S = (1./(d2.shape[0]-1))*d2.T.dot(U2)
				dProj = np.sqrt(d2.shape[0]-1)*S.dot(np.diag(self.evals[[xax-1, yax-2]]**-0.5))
		if not coords:
			f, ax = py.subplots()
			ax.axvline(0, ls='solid', c='k')
			ax.axhline(0, ls='solid', c='k')
			ax.plot(self.U[:,xax-1], self.U[:,yax-1], 'ko', ms=0)
			[ax.text(x,y,s, color=siteCol, fontsize=siteSize, ha='center', va='center') for x,y,s in zip(self.U[:,xax-1], self.U[:,yax-1], self.siteLabs)]
			if descriptors is not None:
				ax.plot(dProj[:,0], dProj[:,1], 'ko', ms=0)
				[ax.text(x,y,s, color=spCol, fontsize=spSize, ha='center', va='center') for x,y,s in zip(dProj[:,0], dProj[:,1], dLabs)]
			ax.set_xlabel('PCoA Axis {!s}'.format(xax))
			ax.set_ylabel('PCoA Axis {!s}'.format(yax))
			py.show()
		else:
			if descriptors is not None:
				return {'Objects': self.U[:,[xax-1, yax-1]], 'Descriptors': dProj}
			else:
				return {'Objects': self.U[:,[xax-1, yax-1]]}

	def shepard(self, xax=1, yax=2):
		coords = self.U[:,[xax-1, yax-1]]
		reducedD = np.zeros((coords.shape[0], coords.shape[0]))
		for i in xrange(coords.shape[0]):
			for j in xrange(coords.shape[0]):
				d = coords[i,:] - coords[j,:]
				reducedD[i, j] = np.sqrt( d.dot(d) )
		reducedD = reducedD[np.tril_indices_from(reducedD, k=-1)]
		originalD = self.y2[np.tril_indices_from(self.y2, k=-1)]
		xmin = np.min(reducedD)
		xmax = np.max(reducedD)
		f, ax = py.subplots()
		ax.plot(reducedD, originalD, 'ko')
		ax.plot([xmin, xmax], [xmin, xmax], 'r--')
		ax.set_xlabel('Distances in Reduced Space')
		ax.set_ylabel('Distances in Original Matrix')
		py.show()




	