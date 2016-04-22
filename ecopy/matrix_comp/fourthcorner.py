import numpy as np
from pandas import DataFrame, factorize
from scipy.stats import chisquare

class corner4(object):
	"""
	Docstring for function ecopy.corner4
	====================
	Conducts a fourth corner analysis for a trait matrix (Q),
		a species matrix (L), and an environmental matrix (R)

	Use
	----
	corner4(R, L, Q, nperm=999, model=1, test='both', p_adjustment=None)

	Returns a pandas.DataFrame containing output summary

	Parameters
	----------
	R: A pandas.DataFrame containing environmental data (site x environmental data)
	L:  A numpy.ndarray or pandas.DataFrame containing species presence/absence
		or abundances at each site (site x species). Only integer values allowed.
	Q: A pandas.DaraFrame containing species traits (species x trait)
	model: Which permutational model to use (1,2,3,4). See online documentation.
	test: Whether the test should be one-tailed ('greater', 'lower') or two-tailed ('both')
	p_adjustment: Multiple comparison adjustment to use for p-values. Can be None,
		'bonferroni', 'holm', or 'fdr'

	Example
	--------
	import ecopy as ep

	traits = ep.load_data('avi_traits')
	env = ep.load_data('avi_env')
	sp = ep.load_data('avi_sp')

	test1 = ep.corner4(env, sp, traits, nperm=99, p_adjustment='fdr')
	print(test1.summary())
	"""
	def __init__(self, R, L, Q, nperm=999, model=1, test='both', p_adjustment=None):
		if not isinstance(R, (DataFrame)):
			msg = 'Matrix R must be a pandas.DataFrame'
			raise ValueError(msg)
		if not isinstance(L, (np.ndarray, DataFrame)):
			msg = 'Matrix L must be a numpy.ndarray or pandas.DataFrame'
			raise ValueError(msg)
		if not isinstance(Q, (DataFrame)):
			msg = 'Matrix Q must be a pandas.DataFrame'
			raise ValueError(msg)
		if model not in [1,2,3,4]:
			msg = 'model must be 1, 2, 3, or 4'
			raise ValueError(msg)
		if isinstance(L, DataFrame):
			L = np.array(L)
		if L.dtype != 'int':
			msg = 'Matrix L must contain only integer values of species abundances (no floats).\nConvert matrix to integers before proceeding'
			raise ValueError(msg)
		if np.any(L < 0):
			msg ='Matrix L cannot have negative values'
			raise ValueError(msg)
		if test not in ['greater', 'lower', 'both']:
			msg = 'test argument must be greater, lower, or both'
			raise ValueError(msg)
		if p_adjustment is not None:
			if p_adjustment not in ['bonferroni', 'holm', 'fdr']:
				msg = 'p_adjustment must be bonferroni, holm, fdr, or None'
				raise ValueError(msg)
		pval = []
		obsStat = []
		statType = []
		compName = []
		for i in range(R.shape[1]):
			for j in range(Q.shape[1]):
				if R.iloc[:,i].dtype not in ['float', 'object']:
					msg = 'Environmental data must be float or object'
					raise ValueError(msg)
				if Q.iloc[:,j].dtype not in ['float', 'object']:
					msg = 'Trait data must be float or object'
					raise ValueError(msg)
				stat_perm = np.zeros(nperm)
				comp = '{0} - {1}'.format(R.iloc[:,i].name, Q.iloc[:,j].name)
				compName.append(comp)
				if R.iloc[:,i].dtype=='float' and Q.iloc[:,i].dtype=='float':
					stat_obs = QuantQuant(R.iloc[:,i], L, Q.iloc[:,j])
					iteration = 0
					while iteration < nperm:
						Lperm = permuteType(L, model)
						stat_perm[iteration] = QuantQuant(R.iloc[:,i], Lperm, Q.iloc[:,j])
						iteration += 1
					obsStat.append(stat_obs)
					statType.append('Pearson r')
					if test is 'greater':
						pval.append(np.mean(stat_perm >= stat_obs))
					if test is 'lower':
						pval.append(np.mean(stat_perm <= stat_obs))
					if test is 'both':
						pval.append(np.mean(np.abs(stat_perm)>=np.abs(stat_obs)))
				if R.iloc[:,i].dtype=='object' and Q.iloc[:,j].dtype=='object':
					R2 = factorize(R.iloc[:,i])[0]
					Q2 = factorize(Q.iloc[:,j])[0]
					n_env = len(R.iloc[:,i].unique())
					n_trait = len(Q.iloc[:,j].unique())
					stat_obs = QualQual(R2, L, Q2, n_env, n_trait)
					iteration = 0 
					while iteration < nperm:
						Lperm = permuteType(L, model)
						stat_perm[iteration] = QualQual(R2, Lperm, Q2, n_env, n_trait)
						iteration += 1
					obsStat.append(stat_obs)
					statType.append('Chi-Squared')
					if test is 'greater':
						pval.append(np.mean(stat_perm >= stat_obs))
					if test is 'lower':
						pval.append(np.mean(stat_perm <= stat_obs))
					if test is 'both':
						pval.append(np.mean(np.abs(stat_perm)>=np.abs(stat_obs)))
				if R.iloc[:,i].dtype=='object' and Q.iloc[:,j].dtype=='float':
					stat_obs = QualQuant(R.iloc[:,i], L, Q.iloc[:,j])
					iteration = 0
					while iteration < nperm:
						Lperm = permuteType(L, model)
						stat_perm[iteration] = QualQuant(R.iloc[:,i], Lperm, Q.iloc[:,j])
						iteration += 1
					obsStat.append(stat_obs)
					statType.append('F')
					if test is 'greater':
						pval.append(np.mean(stat_perm >= stat_obs))
					if test is 'lower':
						pval.append(np.mean(stat_perm <= stat_obs))
					if test is 'both':
						pval.append(np.mean(np.abs(stat_perm)>=np.abs(stat_obs)))
				if R.iloc[:,i].dtype=='float' and Q.iloc[:,j].dtype=='object':
					stat_obs = QuantQual(R.iloc[:,i], L, Q.iloc[:,j])
					iteration = 0
					while iteration < nperm:
						Lperm = permuteType(L, model)
						stat_perm[iteration] = QuantQual(R.iloc[:,i], Lperm, Q.iloc[:,j])
						iteration += 1
					obsStat.append(stat_obs)
					statType.append('F-statistic')
					if test is 'greater':
						pval.append(np.mean(stat_perm >= stat_obs))
					if test is 'lower':
						pval.append(np.mean(stat_perm <= stat_obs))
					if test is 'both':
						pval.append(np.mean(np.abs(stat_perm)>=np.abs(stat_obs)))
		self.results = DataFrame({'Comparison': compName, 'Statistic': statType, 'Observed Stat': np.round(obsStat, 2), 'p-value': np.round(pval, 3)})
		self.results['tail'] = test
		self.nperm = nperm
		self.adj = p_adjustment

	def summary(self):
		if self.adj is None:
			print('\nFour Corner Analysis - {0} Permutations\n'.format(self.nperm))
			return self.results[['Comparison', 'Statistic', 'Observed Stat', 'tail', 'p-value']]
		if self.adj is not None:
			print('\nFour Corner Analysis - {0} Permutations\n{1} Correction\n'.format(self.nperm, self.adj))
			self.results['adjusted p-value'] = p_adjust(self.results['p-value'], self.adj)
			return self.results[['Comparison', 'Statistic', 'Observed Stat', 'tail', 'p-value', 'adjusted p-value']]

def QuantQuant(R, L, Q):
	L2 = L.copy()
	R2 = np.array(R).astype('float').flatten()
	Q2 = np.array(Q).astype('float').flatten()
	Ro = []
	Qo = []
	for i in range(L2.shape[0]):
		for j in range(L2.shape[1]):
			if L2[i,j] != 0:
				Ro.extend([R2[i]]*L2[i,j])
				Qo.extend([Q2[j]]*L2[i,j])
	r = np.corrcoef(np.array(zip(Ro, Qo)), rowvar=0)[0,1]
	return r

def QualQual(R, L, Q, n_env, n_trait):
	L2 = L.copy()
	Ro = []
	Qo = []
	for i in range(L2.shape[0]):
		for j in range(L2.shape[1]):
			if L2[i,j] != 0:
				Ro.extend([R[i]]*L2[i,j])
				Qo.extend([Q[j]]*L2[i,j])
	inflMat = np.array(zip(Ro, Qo))
	crosstab = np.zeros((n_env, n_trait))
	for i in range(n_env):
		for j in range(n_trait):
			crosstab[i,j] = np.argwhere((inflMat[:,0]==i) & (inflMat[:,1]==j)).shape[0]
	return chisquare(crosstab.ravel())[0]

def QualQuant(R, L, Q):
	L2 = L.copy()
	R2 = factorize(R)[0]
	Q2 = Q.copy()
	Ro = []
	Qo = []
	for i in range(L2.shape[0]):
		for j in range(L2.shape[1]):
			if L2[i,j] != 0:
				Ro.extend([R2[i]]*L2[i,j])
				Qo.extend([Q2[j]]*L2[i,j])
	inflMat = DataFrame({'Ro': Ro, 'Qo': Qo})
	btwn = inflMat.groupby('Ro')['Qo'].agg({'mean': lambda x: len(inflMat)*(x.mean() - inflMat['Qo'].mean())**2})
	btwn = btwn['mean'].sum() / (len(btwn)-1)
	within = inflMat.groupby('Ro')['Qo'].agg({'mean': lambda x: np.sum((x - x.mean())**2)})
	within = within['mean'].sum() / (len(inflMat)-len(within))
	return btwn / within

def QuantQual(R, L, Q):
	L2 = L.copy()
	R2 = R.copy()
	Q2 = factorize(Q)[0]
	Ro = []
	Qo = []
	for i in range(L2.shape[0]):
		for j in range(L2.shape[1]):
			if L2[i,j] != 0:
				Ro.extend([R2[i]]*L2[i,j])
				Qo.extend([Q2[j]]*L2[i,j])
	inflMat = DataFrame({'Ro': Ro, 'Qo': Qo})
	btwn = inflMat.groupby('Qo')['Ro'].agg({'mean': lambda x: len(inflMat)*(x.mean() - inflMat['Qo'].mean())**2})
	btwn = btwn['mean'].sum() / (len(btwn)-1)
	within = inflMat.groupby('Qo')['Ro'].agg({'mean': lambda x: np.sum((x - x.mean())**2)})
	within = within['mean'].sum() / (len(inflMat)-len(within))
	return btwn / within

def permuteType(L, model):
	if model==1:
		Lperm = np.apply_along_axis(lambda x: np.random.permutation(x), 0, L)
	if model==2:
		idx = np.arange(L.shape[0])
		Lperm = L[np.random.permutation(idx),:]
	if model==3:
		Lperm = np.apply_along_axis(lambda x: np.random.permutation(x), 1, L)
	if model==4:
		idx = np.arange(L.shape[1])
		Lperm = L[:,np.random.permutation(idx)]
	return Lperm

def p_adjust(p, method):
	if method == 'bonferroni':
		return np.minimum(p*len(p), 1)
	if method == 'holm':
		temp = DataFrame({'p': p})
		temp.sort(columns='p', inplace=True)
		temp['newID'] = range(1, len(temp)+1)
		temp['p_adj'] = np.minimum(temp['p'] * (1 + len(temp) - temp['newID']), 1)
		temp.sort(inplace=True)
		return temp['p_adj']
	if method == 'fdr':
		temp = DataFrame({'p': p})
		temp.sort(columns='p', inplace=True, ascending=False)
		temp['newID'] = range(1, len(temp)+1)
		temp['p_adj'] = np.minimum(1, len(temp)/temp['newID'] * temp['p'])
		temp.sort(inplace=True)
		return np.round(temp['p_adj'], 3)
