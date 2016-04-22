import numpy as np
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import pearsonr
from ..regression import isotonic
from ..ordination import pca, pcoa

class MDS(object):
	'''
	Docstring for function ecopy.pca
	====================
	Conducts multidimensional scaling  (MDS). User
		supplies a square-symmetric distance matrix

	Use
	----
	MDS(distmat, siteNames=None, naxes=2, transform='monotone', ntry=20, tolerance=1E-4, maxiter=3000, init=None)

	Returns an object of class MDS

	Parameters
	----------
	distmat: a square-symmetric distance matrix. Can be either a numpy.ndarray or pandas.DataFrame
	siteNames: A list of names for each object. If none, takes on integer values or the index of the pandas.DataFrame
	naxes: Number of ordination axes. Default = 2
	transform: Which transformation should be used during scaling.
		'absolute': Conducts absolute MDS. Distances between points in ordination space
			should be as close as possible to observed distances.
		'ratio': Ordination distances are proportional to observed distances.
		'linear': Ordination distances are a linear function of observed distances. Uses the 
			technique of Heiser (1991) to avoid negative ordination distances.
		'monotone':  Constrains ordination distances simply to be ranked the same as 
			observed distance. Typically referred to as non-metric multidimensional 
			scaling. Uses isotonic regression from scikit-learn. Default value.
	ntry: Number of random starts used to avoid local minima. The returned solution is
		the one with the lowest final stress.
	tolerance: Minimum step size causing a break in the minimization of stress. Default = 1E-4.
	maxiter: Maximum number of iterations to attempt before breaking if no solution is found.
	init: Initial positions for the first random start. If none, the initial position of the first try is taken
		as the site locations from classical scaling, Principle Coordinates Analysis.

	See online documentation for more details.

	Attributes
	---------
	scores: Final object scores on the ordination axes.
	stress: Final stress.
	obs: The observed distance matrix.
	transform: Which transformation was used.

	Methods
	--------
	biplot(xax=1, yax=2, siteNames=True, descriptors=None, descripNames=None, coords=False):
		A biplot of ordination axes.
			xax: Which ordination axis should be on the x-axis
			yax: Which ordination axis should be on the y-axis
			siteNames: Whether sites should be plotted as text (True) or points (False)
			descriptors: A numpy.ndarray or pandas.DataFrame of predictors used to construct
				the distance matrix (i.e. species). Species scores are the weighted average
				along each axis. 
			descripNames: A list of names for the descriptors. If None, uses column names
				 from the pandas.DataFrame.  Otherwise, uses integers.
			coords: If True, no plot is produced, but a dictionary of Object and Descriptor scores
				returned for custom plotting
	shepard(): Produces a shepard plot of observed vs. fitted distances.
	correlations(): Returns the correlation between observed and fitted distances
		for each site.
	correlationPlots(site=None): Produces a plot of observed vs. fitted distances for a given
		site. If site is None, all sites are plotted in a single graph.

	Example
	--------
	import ecopy as ep
	dunes = ep.load_data('dune')
	dunes_T = ep.transform(dunes, 'wisconsin')
	dunes_D = ep.distance(dunes_T, 'bray')
	dunesMDS = ep.MDS(dunes_D, transform='monotone')
	dunesMDS.biplot()
	dunesMDS.biplot(descriptors=dunes_T)
	'''
	def __init__(self, distmat, siteNames=None, naxes=2, transform='monotone', 
			ntry=20, tolerance=1E-4, maxiter=3000, init=None):
		if transform not in ['monotone', 'absolute', 'linear', 'ratio']:
			msg = 'transform must be one of monotone, absolute, linear, ratio'
			raise ValueError(msg)
		if not isinstance(distmat, (np.ndarray, DataFrame)):
			msg = 'distmat must be a square symmetric numpy.ndarray or pandas.DataFrame'
			raise ValueError(msg)
		if isinstance(distmat, DataFrame):
			if (distmat.dtypes == 'object').any():
				msg = 'DataFrame can only contain numeric values'
				raise ValueError(msg)
		distmat = np.array(distmat)
		if distmat.any() < 0:
			msg ='Distance matrix cannot contain negative values'
			raise ValueError(msg)
		if distmat.shape[0] != distmat.shape[1]:
			msg = 'distmat must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if not np.allclose(distmat.T, distmat):
			msg ='distmat must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if siteNames is not None:
			self.siteLabs = siteNames
		else:
			self.siteLabs = ['Site ' + str(x) for x in range(1, distmat.shape[0]+1)]
		weights = np.ones(distmat.shape)
		weights[np.isnan(distmat)] = 0
		coordStore = None
		stressStore = 1
		if init is None:
			init2 = pcoa(distmat).U[:,:naxes]
		else:
			init2 = init
		if transform is 'absolute':
			Vp = VTrans(weights)
			for i in range(ntry):
				if i ==0:
					Z = init2
				else:
					Z = np.random.rand(distmat.shape[0]*naxes).reshape(distmat.shape[0], naxes)
				stress1 = 1
				k = 0
				e = 1E4
				while k < maxiter and e > tolerance:
					stress2, Xu = absMDS(distmat, Z, weights, Vp)
					e = stress1 - stress2
					Z = Xu
					stress1 = stress2
					k += 1
				print('Finished at iteration {0}. Stress = {1}'.format(k, stress2))	
				if stress2 < stressStore:
					stressStore = stress2
					coordStore = Z
				self.parameters = None
		if transform is 'ratio':
			bStore = None
			Vp = VTrans(weights)
			for i in range(ntry):
				if i ==0:
					Z = init2
				else:
					Z = np.random.rand(distmat.shape[0]*naxes).reshape(distmat.shape[0], naxes)
				stress1 = 1
				k = 0
				e = 1E4
				b = np.random.rand(1)
				while k < maxiter and e > tolerance:
					stress2, Xu, b = ratioMDS(distmat, b, Z, weights, Vp)
					e = stress1 - stress2
					Z = Xu
					stress1 = stress2
					k += 1
				print('Finished at iteration {0}. Stress = {1}'.format(k, stress2))
				if stress2 < stressStore:
					stressStore = stress2
					coordStore = Z
					bStore = b
				self.parameters = {'b': bStore}
		if transform is 'linear':
			aStore = None
			bStore = None
			for i in range(ntry):
				if i ==0:
					Z = init2
				else:
					Z = np.random.rand(distmat.shape[0]*naxes).reshape(distmat.shape[0], naxes)
				stress1 = 1
				k = 0
				e = 1E4
				a = np.random.rand(1)
				b = np.random.rand(1)
				while k < maxiter and e > tolerance:
					stress2, Xu, a, b, = linMDS(distmat, a, b, Z)
					e = stress1 - stress2
					Z = pca(Xu).scores
					stress1 = stress2
					k += 1
				print('Finished at iteration {0}. Stress = {1}'.format(k, stress2))
				if stress2 < stressStore:
					stressStore = stress2
					coordStore = Z
					aStore = a
					bStore = b
				self.parameters = {'a': aStore, 'b': bStore}
		if transform is 'monotone':
			Vp = VTrans(weights)
			for i in range(ntry):
				if i ==0:
					Z = init2
				else:
					Z = np.random.rand(distmat.shape[0]*naxes).reshape(distmat.shape[0], naxes)
				stress1 = 1
				k = 0
				e = 1E4
				while k < maxiter and e > tolerance:
					stress2, Xu = nMDS(distmat, Z, weights, Vp)
					e = stress1 - stress2
					Z = pca(Xu).scores
					stress1 = stress2
					k += 1
				print('Finished at iteration {0}. Stress = {1}'.format(k, stress2))
				if stress2 < stressStore:
					stressStore = stress2
					coordStore = Z
		self.scores = np.apply_along_axis(lambda x: x - np.nanmean(x), 0, coordStore)
		self.stress = np.min(stressStore)
		self.scores = scoreTrans(self.scores, distmat)
		self.scores = pca(self.scores).scores
		self.obs = distmat
		self.transform = transform
		print('Final Stress = {0}'.format(self.stress))

	def biplot(self, xax=1, yax=2, siteNames=True, coords=False, descriptors=None, descripNames=None, spCol='r', siteCol='k', spSize=12, siteSize=12):
		if descriptors is not None:
			wts = np.diag(descriptors.sum(axis=0)**-1)
			dScores = wts.dot(descriptors.T).dot(self.scores)
			if descripNames is None:
				if isinstance(descriptors, DataFrame):
					descripNames = descriptors.columns
				else:
					descripNames = ['Species {0}'.format(i) for i in range(descriptors.shape[1])]
		if not coords:
			f, ax = plt.subplots()
			ax.axvline(0, ls='solid', c='k')
			ax.axhline(0, ls='solid', c='k')
			if siteNames:
				ax.plot(self.scores[:,xax-1], self.scores[:,yax-1], 'ko', ms=0)
				[ax.text(x,y,s, color=siteCol, fontsize=siteSize, ha='center', va='center') for x,y,s in zip(self.scores[:,xax-1], self.scores[:,yax-1], self.siteLabs)]
			else:
				ax.plot(self.scores[:,xax-1], self.scores[:,yax-1], 'ko', ms=8)
			if descriptors is not None:
				[ax.text(x,y,s, color=spCol, fontsize=spSize, ha='center', va='center') for x,y,s in zip(dScores[:,0], dScores[:,1], descripNames)]
			ax.set_xlabel('nMDS Axis {!s}'.format(xax))
			ax.set_ylabel('nMDS Axis {!s}'.format(yax))
			plt.show()
		else:
			if descriptors is None:
				coordDict = {'Objects': self.scores}
			else:
				coordDict = {'Objects': self.scores, 'Descriptors': dScores}
			return coordDict

	def shepard(self):
		dHats = np.tril(eucD(self.scores)).ravel()
		dObs = np.tril(self.obs).ravel()
		f, ax = plt.subplots()
		ax.plot(dObs, dHats, 'ro')
		if self.transform is 'linear':
			X = np.column_stack((np.ones(len(dObs)), dObs))
			Y = dHats
			B = np.linalg.solve(X.T.dot(X), X.T.dot(Y))
			linPreds = B[0] + B[1]*dObs
			R2 = 1 - ((dHats - linPreds)**2).sum() / ((dHats-dObs)**2).sum()
			o = dObs.argsort()
			ax.plot(dObs[o], linPreds[o], c='b', lw=2)
			ax.text(0.05, 0.85, 'R2 = {:.3}'.format(R2), ha='left', va='center', transform=ax.transAxes)
		if self.transform is 'monotone':
			monoPreds = isotonic(dHats, dObs)
			R2 = 1 - ((dHats - monoPreds.prediction)**2).sum() / ((dHats-dObs)**2).sum()
			o = dObs.argsort()
			ax.step(dObs[o], monoPreds.prediction[o], c='b', lw=2)
			ax.text(0.05, 0.85, 'R2 = {:.3}'.format(R2), ha='left', va='center', transform=ax.transAxes)
		ax.text(0.05, 0.9, 'Stress = {:.3}'.format(self.stress), ha='left', va='center', transform = ax.transAxes)
		plt.show()

	def correlations(self):
		dFit = eucD(self.scores)
		dObs = self.obs
		corrs = []
		for i in range(len(dObs)):
			corrs.append(pearsonr(dObs[0,:], dFit[0,:])[0])
		corrs = Series(corrs)
		return corrs

	def correlationPlots(self, site=None):
		f, ax = plt.subplots()
		dFit = eucD(self.scores)
		dObs = self.obs
		if site is None:
			for i in range(len(self.obs)):
				ax.plot(dObs[i,:], dFit[i,:], 'o', alpha=0.5, c=cm.Set1(float(i)/len(dObs), 1))
		else:
			ax.plot(dObs[site,:], dFit[site,:], 'ko', ms=0)
			[ax.text(x,y,str(s), ha='center', va='center') for x,y,s in zip(dObs[site,:], dFit[site,:], range(len(dFit[site,:])))]
			ax.set_title('Site {0}'.format(site))
		ax.set_ylabel("Fitted Distance to Site")
		ax.set_xlabel("Observed Distance to Site")
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		plt.show()
		

def VTrans(weights):
	V = - weights
	np.fill_diagonal(V, 'nan')
	np.fill_diagonal(V, np.nansum(-1*V, axis=1))
	return np.linalg.pinv(V)

def Bcalc(weights, distmat, dZ):
	ratio = weights*distmat/dZ
	ratio[dZ==1E-5] = 0
	bZ = -ratio
	np.fill_diagonal(bZ, ratio.sum(axis=1))
	return bZ

def eucD(X):
	n = len(X)
	one = np.ones((n, 1))
	D = np.sum(X**2, 1).reshape(n, 1).dot(one.T) + one.dot(np.sum(X**2, 1).reshape(1, n)) - 2*X.dot(X.T)
	D = np.sqrt(D)
	return D

def absMDS(distmat, Z, weights, Vp):
	dZ = eucD(Z)
	dZ[dZ==0] = 1E-5
	bZ = Bcalc(weights, distmat, dZ)
	Xu = Vp.dot(bZ).dot(Z)
	dXu = eucD(Xu)
	stress = np.sqrt(np.tril(weights*(distmat-dXu)**2).sum() / np.tril(dXu**2).sum())
	return stress, Xu

def ratioMDS(distmat, b, Z, weights, Vp):
	dHat = distmat*b
	dZ = eucD(Z)
	dZ[dZ==0] = 1E-5
	bZ = Bcalc(weights, dHat, dZ)
	Xu = Vp.dot(bZ).dot(Z)
	dXu = eucD(Xu)
	stress = np.sqrt(np.tril(weights*(dHat-dXu)**2).sum() / np.tril(dXu**2).sum())
	b = np.tril(weights*distmat*dXu).sum() / np.tril(weights*distmat**2).sum()
	return stress, Xu, b

def linMDS(distmat, a, b, Z):
	dHat = a + b*distmat
	dZ = eucD(Z)
	dZ[dZ==0] = 1E-5
	weights = np.ones(distmat.shape)
	weights[np.isnan(distmat)] = 0 
	ix = dZ<0
	weights[ix] = weights[ix]*(dZ[ix] + np.abs(dZ[ix]))/dZ[ix]
	Vp = VTrans(weights)
	bZ = Bcalc(weights, dHat, dZ)
	bZ[ix] = 0
	Xu = Vp.dot(bZ).dot(Z)
	dXu = eucD(Xu)
	stress = np.sqrt(np.tril(weights*(dHat-dXu)**2).sum() / np.tril(dXu**2).sum())
	a = (np.tril(weights*dXu).sum() - b*np.tril(weights*distmat).sum()) / np.tril(weights).sum()
	b = (np.tril(weights*distmat*dXu).sum() - a*np.tril(weights*distmat).sum()) / np.tril(weights*distmat**2).sum()
	return stress, Xu, a, b

def nMDS(distmat, Z, weights, Vp):
	dZ = eucD(Z)
	dZ[dZ==0] = 1E-5
	diss_f = distmat.ravel()
	dhat_f = dZ.ravel()
	dhat = isotonic(dhat_f, diss_f)
	dhat = dhat.prediction.reshape(distmat.shape)
	stress = np.sqrt(np.tril((weights*(dZ - dhat)**2)).sum() / np.tril(weights*dZ**2).sum())
	bZ = Bcalc(weights, dhat, dZ)
	Xu = Vp.dot(bZ).dot(Z)
	return stress, Xu

def scoreTrans(x, distmat):
	distMax = eucD(x).max()
	obsMax = x.max()
	scale = obsMax / distMax
	return x*scale
