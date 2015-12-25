import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt

class anosim(object):
	'''
	Docstring for function ecopy.anosim
	====================
	Conducts analysis of similarity (ANOSIM) on a distance matrix given
		one or two factors (groups)

	Use
	----
	anosim(dist, factor1, factor2=None, nested=False, nperm=999)

	Returns an object of class anosim

	Parameters
	----------
	dist:  A square-symmetric distance matrix
	factor1: The first grouping factor
	factor2: The second grouping factor
 	nested: Whether or not the first factor is nested within the second.
 		If true, levels of factor1 are permuted only within the nesting groups
 	nperm: Number of permutations for the test.

	Attributes (see online documentation for descriptions)
	---------
	r_perm1: Permuted R-statistics of the null distribution for factor1
	r_perm2: Permuted R-statistics of the null distribution for factor2
	R_obs1: The observed R-statistic of factor1
	R_obs2: The observed R-statistic of factor2
	pval: List of pvals for factor1 and factor2
	perm: Number of permutations

	Methods
	--------
	summary(): provides a summary of test results
	plot(): plots histograms of random vs. observed
		R-statistics

	Example
	--------
	import ecopy as ep

	data1 = ep.load_data('dune')
	data2 = ep.load_data('dune_env')
	duneDist = ep.distance(data1, 'bray')
	group1 = data2['Management']
	group2map = {'SF': 'A', 'BF': 'A', 'HF': 'B', 'NM': 'B'}
	group2 = group1.map(group2map)
	t1 = ep.anosim(duneDist, group1, group2, nested=True, nperm=9999)
	print(t1.summary())
	t1.plot()
	'''
	def __init__(self, dist, factor1, factor2=None, nested=False, nperm=999):
		if isinstance(dist, DataFrame):
			dist = np.array(dist)
		if dist.shape[0] != dist.shape[1]:
			msg = 'Matrix dist must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if not np.allclose(dist.T, dist):
			msg = 'Matrix dist must be a square, symmetric distance matrix'
			raise ValueError(msg)
		if np.any(dist < 0):
			msg = 'Distance matrix cannot have negative values'
			raise ValueError(msg)
		self.r_perm1 = np.empty(nperm)
		self.r_perm2 = np.empty(nperm)
		self.R_obs1 = None
		self.R_obs2 = None
		if factor2 is None:
			g1 = np.array(factor1)
			self.R_obs1 = oneWayANOSIM(dist, g1)
			for i in range(nperm):
				groupRand = np.random.choice(g1, len(g1), replace=False)
				self.r_perm1[i] = oneWayANOSIM(dist, groupRand)
			self.p_val = np.mean(self.r_perm1 > self.R_obs1)
		if factor2 is not None and not nested:
			g1 = np.array(factor1)
			self.R_obs1 = oneWayANOSIM(dist, g1)
			for i in range(nperm):
				groupRand = np.random.choice(g1, len(g1), replace=False)
				self.r_perm1[i] = oneWayANOSIM(dist, groupRand)
			g2 = np.array(factor2)
			self.R_obs2 = oneWayANOSIM(dist, g2)
			for i in range(nperm):
				groupRand = np.random.choice(g2, len(g2), replace=False)
				self.r_perm2[i] = oneWayANOSIM(dist, groupRand)
			self.p_val = [np.mean(self.r_perm1 > self.R_obs1), np.mean(self.r_perm2 > self.R_obs2)]
		if factor2 is not None and nested:
			g1 = np.array(factor1)
			g2 = np.array(factor2)
			comb = np.array(zip(g1, g2), dtype=[('group1', 'S10'), ('group2', 'S10')])
			sortIX = comb.argsort(order='group2')
			comb = comb[sortIX]
			dist1 = dist[sortIX,:][:,sortIX]
			gpR = []
			for i in np.unique(comb['group2']):
				withinMat = dist1[comb['group2']==i,:][:,comb['group2']==i]
				gpR.append(oneWayANOSIM(withinMat, comb['group1'][comb['group2']==i]))
			self.R_obs1 = np.mean(gpR)
			comb2 = comb.copy()
			for i in range(nperm):
				permG = []
				for j in np.unique(comb2['group2']):
					gPerm = np.random.choice(comb2['group1'][comb2['group2']==j], np.sum(comb2['group2']==j), replace=False)
					permG.extend(gPerm)
				comb2['group1'] = permG
				gpR = []
				for k in np.unique(comb2['group2']):
					withinMat = dist1[comb['group2']==k,:][:,comb['group2']==k]
					gpR.append(oneWayANOSIM(withinMat, comb2['group1'][comb2['group2']==k]))
				self.r_perm1[i] = np.mean(gpR)
			dist2 = dist1.copy()
			li = np.tril_indices(dist1.shape[0])
			dist2[li] = np.nan
			rankMat = dist2.flatten().argsort().argsort().reshape(dist2.shape)
			uniqueSites = np.unique(comb)[np.unique(comb).argsort(order='group2')]['group1']
			collapseMat =np.zeros((len(uniqueSites), len(uniqueSites)))
			for i in range(len(uniqueSites)):
				for j in range(len(uniqueSites)):
					collapseMat[i,j] = np.mean(rankMat[comb['group1']==uniqueSites[i],:][:,comb['group1']==uniqueSites[j]])
			np.fill_diagonal(collapseMat, 0)
			collapseGroup = np.unique(comb)['group2']
			self.R_obs2 = oneWayANOSIM(collapseMat, collapseGroup)
			for i in range(nperm):
				groupRand = np.random.choice(collapseGroup, len(collapseGroup), replace=False)
				self.r_perm2[i] = oneWayANOSIM(collapseMat, groupRand)
			self.p_val = [np.mean(self.r_perm1 > self.R_obs1), np.mean(self.r_perm2 > self.R_obs2)]
		self.perm = nperm


	def summary(self):
		if self.R_obs2 is None:
			summ1 = '\nANOSIM\nObserved R = {0:.3}\np-value = {1:.3}\n{2} permutations'.format(self.R_obs1, self.p_val, self.perm)
			return summ1
		else:
			summ1 = '\nANOSIM: Factor 1\nObserved R = {0:.3}\np-value = {1:.3}\n{2} permutations'.format(self.R_obs1, self.p_val[0], self.perm)
			summ2 = '\nANOSIM: Factor 2\nObserved R = {0:.3}\np-value = {1:.3}\n{2} permutations'.format(self.R_obs2, self.p_val[1], self.perm)
			return summ1 + '\n' + summ2
	def plot(self):
		if self.R_obs2 is None:
			f, ax = plt.subplots()
			ax.hist(self.r_perm1, 50, normed=1, color='blue', alpha=0.5, histtype='stepfilled', linewidth=1)
			ax.axvline(self.R_obs1, linewidth=2, linestyle='dashed', color='red', label='Observed R')
			ax.set_ylabel("Density")
			ax.set_xlabel("R-statistic")
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_ticks_position('bottom')
			ax.legend(loc=1)
			plt.show()
		else:
			f, ax = plt.subplots(2, 1, figsize=(6.5, 8.5))
			ax[0].hist(self.r_perm1, 50, normed=1, color='blue', alpha=0.5, histtype='stepfilled', linewidth=1)
			ax[0].axvline(self.R_obs1, linewidth=2, linestyle='dashed', color='red', label='Observed R')
			ax[0].set_ylabel("Density")
			ax[0].set_xlabel("R-statistic")
			ax[0].spines['top'].set_visible(False)
			ax[0].spines['right'].set_visible(False)
			ax[0].yaxis.set_ticks_position('left')
			ax[0].xaxis.set_ticks_position('bottom')
			ax[0].legend(loc=1)
			ax[0].set_title('Factor 1')

			ax[1].hist(self.r_perm2, 50, normed=1, color='red', alpha=0.5, histtype='stepfilled', linewidth=1)
			ax[1].axvline(self.R_obs2, linewidth=2, linestyle='dashed', color='red', label='Observed R')
			ax[1].set_ylabel("Density")
			ax[1].set_xlabel("R-statistic")
			ax[1].spines['top'].set_visible(False)
			ax[1].spines['right'].set_visible(False)
			ax[1].yaxis.set_ticks_position('left')
			ax[1].xaxis.set_ticks_position('bottom')
			ax[1].legend(loc=1)
			ax[1].set_title('Factor 2')
			plt.show()

def oneWayANOSIM(x, group):
	sortIX = group.argsort()
	sorted_G = group[sortIX]
	sorted_Dist = x[sortIX,:][:,sortIX]
	groupMat = np.array([sorted_G == i for i in sorted_G])
	ui = np.triu_indices(groupMat.shape[0], k=1)
	distU = sorted_Dist[ui]
	distU = distU.argsort().argsort()
	groupU = groupMat[ui]
	r_w = np.mean(distU[groupU==True])
	r_b = np.mean(distU[groupU==False])
	n = x.shape[0]
	denom = n*(n-1)/4.
	return (r_b - r_w) / denom	
