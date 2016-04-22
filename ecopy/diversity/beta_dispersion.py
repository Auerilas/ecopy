import numpy as np
from pandas import DataFrame, Series
from scipy.stats import f
from ..base_funcs import spatial_median

def beta_dispersion(X, groups, test='anova', scores=False, center='median', n_iter=99):
	'''
	Docstring for function ecopy.beta_dispersion
	========================
	Calculates beta dispersion among groups for a given distance matrix.

	Use
	----
	beta_dispersion(x, method='shannon', breakNA=True)

	Returns an ANOVA table for beta dispersion. If scores=True,
	returns a numpy.ndarray of calculated z-distances (see online
	documentation)

	Parameters
	----------
	X:  numpy.ndarray or pandas.DataFrame containing a square
		symmetric distance matrix
	groups: a list, pandas.Series, or pandas.DataFrame containing the
		group for each observation
	test: either 'anova' or 'permute' for calculating sigficiance
	scores: Boolean, whether or not the calculated z-distances should
		be returned for plotting
	center: Either 'median' or 'centroid', identifies which mid-point
		should be used to calculate z-distances
	n_iter: integer, number of permutations

	Example
	--------
	import ecopy as ep
	varespec = ep.load_data('varespec')
	dist = ep.distance(varespec, method='bray')

	groups = ['grazed']*16 + ['ungrazed']*8
	print(ep.beta_dispersion(dist, groups, test='permute', center='median', scores=False))
	'''
	if not isinstance(X, (np.ndarray, DataFrame)):
		msg = 'X must be a numpy.ndarray or pandas.DataFrame'
		raise ValueError(msg)
	if not isinstance(groups, (list, Series, DataFrame)):
		msg = 'groups must be a list, pandas.Series, or pandas.DataFrame'
		raise ValueError(msg)
	Y = np.array(X)
	groups = list(groups)
	if np.isnan(Y).any():
		msg = 'Distance matrix contains null values'
		raise ValueError(msg)
	if Y.any() < 0:
		msg ='Distance matrix cannot contain negative values'
		raise ValueError(msg)
	if Y.shape[0] != Y.shape[1]:
		msg = 'Distance matrix must be square'
		raise ValueError(msg)
	if not np.allclose(Y.T, Y):
		msg ='Distance matrix must be symmetric'
		raise ValueError(msg)
	if len(groups) != Y.shape[0]:
		msg = 'Number of groups does not equal the rows of the distance matrix'
		raise ValueError(msg)
	if center not in ['median', 'centroid']:
		msg = 'center must be either median or centroid'
		raise ValueError(msg)
	if test not in ['anova', 'permute']:
		msg = 'test must be either anova or permute'
		raise ValueError(msg)
	A = -0.5*Y**2
	r_means = A.mean(axis=1)
	c_means = A.mean(axis=0)
	o_mean = A.mean()
	G = A - r_means[:,np.newaxis] - c_means[np.newaxis,:] + o_mean
	evals, U = np.linalg.eig(G)
	idx = np.argsort(evals)[::-1]
	evals = evals[idx]
	U = U[:,idx]
	Usc = U.dot(np.diag(np.abs(evals)**0.5))
	if test=='anova':
		obs = z_calc(Usc, groups, evals, center, n_iter)
		print('{0:<10} {1:^5} {2:^8} {3:^8} {4:^8} {5:^8}'.format(' ','df', 'SS', 'MS', 'F', 'P(>F)'))
		print('{0:<10} {1:^5} {2:^8.4f} {3:^8.4f} {4:^8.4f} {5:^8.4f}'.format('Groups', int(obs[4]), obs[0], obs[2], obs[6], 1-f.cdf(obs[6], obs[4], obs[5])))
		print('{0:<10} {1:^5} {2:^8.4f} {3:^8.4f}'.format('Residuals', int(obs[5]), obs[1], obs[3]))
	if test=='permute':
		obs = z_calc(Usc, groups, evals, center, permute=True, iterations=n_iter)
		print('{0:<10} {1:^5} {2:^8} {3:^8} {4:^8} {5:^8}'.format(' ','df', 'SS', 'MS', 'F', 'P(>F)'))
		print('{0:<10} {1:^5} {2:^8.4f} {3:^8.4f} {4:^8.4f} {5:^8.4f}'.format('Groups', int(obs[4]), obs[0], obs[2], obs[6][0], np.mean(obs[6] > obs[6][0])))
		print('{0:<10} {1:^5} {2:^8.4f} {3:^8.4f}'.format('Residuals', int(obs[5]), obs[1], obs[3]))
	if scores:
		return obs[7]


def z_calc(U, groups, evals, center, iterations, permute=False):
	groupIDs = list(set(groups))
	n_groups = len(groupIDs)
	n = len(groups)
	means_groups = np.empty(n_groups)
	n_wg = np.empty(n_groups)
	z_final = np.empty(n)
	res = np.empty(U.shape)
	pos = evals>=0
	for i in range(n_groups):
		idx = np.array(groups)==groupIDs[i]
		pos_Z = U[:,pos][idx,:]
		if center=='median':
			pos_cent = spatial_median(pos_Z)
		if center=='centroid':
			pos_cent = pos_Z.mean(axis=0)
		pos_euclidean = np.apply_along_axis(lambda t: np.sum((t-pos_cent)**2), 1, pos_Z)
		res[np.ix_(idx, pos)] = pos_Z - pos_cent
		neg_euclidean = 0
		neg_res = 0
		neg_cent = 0
		if np.sum(~pos) > 0:
			neg_Z = U[:,~pos][idx,:]
			if center=='median':
				neg_cent = spatial_median(neg_Z)
			if center=='centroid':
				neg_cent = neg_Z.mean(axis=0)
			neg_euclidean = np.apply_along_axis(lambda t: np.sum((t-neg_cent)**2), 1, neg_Z)
			neg_res = neg_Z - neg_cent[np.newaxis,:]
			res[np.ix_(idx, ~pos)] = neg_Z - neg_cent
		temp_r = np.sqrt(np.abs(pos_euclidean-neg_euclidean))
		means_groups[i] = np.mean(temp_r)
		n_wg[i] = len(temp_r)
		z_final[idx] = temp_r
	if not permute:
		ovMean = np.sum(n_wg*means_groups)/np.sum(n_wg)
		SS_bg = np.sum(n_wg*(means_groups-ovMean)**2)
		SS_tot = np.sum((z_final - ovMean)**2)
		SS_wg = SS_tot - SS_bg
		df_bg = n_groups - 1.
		df_wg = float(n-n_groups)
		MS_bg = SS_bg / (df_bg)
		MS_wg = SS_wg / (df_wg)
		F = MS_bg / MS_wg
		return SS_bg, SS_wg, MS_bg, MS_wg, df_bg, df_wg, F, z_final
	if permute:
		ovMean = np.sum(n_wg*means_groups)/np.sum(n_wg)
		SS_bg = np.sum(n_wg*(means_groups-ovMean)**2)
		SS_tot = np.sum((z_final - ovMean)**2)
		SS_wg = SS_tot - SS_bg
		df_bg = n_groups - 1.
		df_wg = float(n-n_groups)
		MS_bg = SS_bg / (df_bg)
		MS_wg = SS_wg / (df_wg)
		F = np.empty(iterations)
		F[0] = MS_bg / MS_wg
		k = 1
		while k < iterations:
			index1 = np.arange(res.shape[0])
			res_perm_pos = res[np.ix_(np.random.choice(index1, len(index1), replace=False),pos)]
			res_perm_pos = res_perm_pos + pos_cent[np.newaxis,:]
			res_perm_neg = res[np.ix_(np.random.choice(index1, len(index1), replace=False),~pos)]
			res_perm_neg = res_perm_neg + neg_cent[np.newaxis,:]
			U_perm = np.concatenate((res_perm_pos, res_perm_neg), axis=1)
			for i in range(n_groups):
				idx = np.array(groups)==groupIDs[i]
				pos_Z = U_perm[np.ix_(idx, pos)]
				pos_euclidean = np.apply_along_axis(lambda t: np.sum((t-pos_cent)**2), 1, pos_Z)
				neg_euclidean = 0
				if np.sum(~pos) > 0:
					neg_Z = U_perm[np.ix_(idx,~pos)]
					neg_euclidean = np.apply_along_axis(lambda t: np.sum((t-neg_cent)**2), 1, neg_Z)
				temp_r = np.sqrt(np.abs(pos_euclidean-neg_euclidean))
				means_groups[i] = np.mean(temp_r)
				n_wg[i] = len(temp_r)
				z_final[idx] = temp_r
			ovMean_perm = np.sum(n_wg*means_groups)/np.sum(n_wg)
			SS_bg_perm = np.sum(n_wg*(means_groups-ovMean_perm)**2)
			SS_tot_perm = np.sum((z_final - ovMean_perm)**2)
			SS_wg_perm = SS_tot_perm - SS_bg_perm
			df_bg = n_groups - 1.
			df_wg = float(n-n_groups)
			MS_bg_perm = SS_bg_perm / (df_bg)
			MS_wg_perm = SS_wg_perm / (df_wg)
			F[k] = MS_bg_perm / MS_wg_perm
			k+=1
		return SS_bg, SS_wg, MS_bg, MS_wg, df_bg, df_wg, F, z_final









