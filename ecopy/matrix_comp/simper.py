import numpy as np
from pandas import DataFrame

def simper(data, factor, spNames=None):
	'''
	Docstring for function ecopy.simper
	====================
	Conducts a SIMPER (percentage similarity) analysis for a 
		site x species matrix given a grouping factor

	Use
	----
	simper(data, factor, spNames=None)

	Returns a pandas.DataFrame containing output

	Parameters
	----------
	data: A numpy.ndarray or pandas.DataFrame site x species matrix
	factor: A numpy.ndarray, pandas.DataFrame, pandas.Series, or list
		containing the factor groupings
	spNames: A list of species names. If data is a pandas.DataFrame,
		spNames are automatically inferred from column names, else
		assigned an integer value.

	Example
	--------
	import ecopy as ep

	data1 = ep.load_data('dune')
	data2 = ep.load_data('dune_env')
	group1 = data2['Management']
	fd = ep.simper(np.array(data1), group1, spNames=data1.columns)
	print(fd.ix['BF-NM'])
	'''
	if not isinstance(data, (DataFrame, np.ndarray)):
		msg = 'datamust be either numpy array or dataframe'
		raise ValueError(msg)
	if isinstance(data, DataFrame):
		if (data.dtypes == 'object').any():
			msg = 'DataFrame can only contain numeric values'
			raise ValueError(msg)
		X = np.array(data).astype('float')
		spNames = data.columns
	else:
		X = data.astype('float')
	if np.min(np.sum(X, axis=1))==0:
		msg = 'One row is entirely zeros, distance calculations will be meaningless'
		raise ValueError(msg)
	if (X < 0).any():
		msg = 'Matrix contains negative values'
		raise ValueError(msg)
	if spNames is None:
		spNames = np.arange(X.shape[1])
	s1 = np.array(spNames)
	g1 = np.array(factor)
	if len(g1) != X.shape[0]:
		msg = 'Factor length must equal number of rows in matrix'
		raise ValueError(msg)
	if len(s1) != X.shape[1]:
		msg = 'Species names must equal number of columns in matrix'
		raise ValueError(msg)
	groupIDs = np.unique(g1)
	print('\nComparison indices:')
	i = 0
	while i < len(groupIDs)-1:
		t1 = X[g1==groupIDs[i],:]
		j = i+1
		while j < len(groupIDs):
			t2 = X[g1==groupIDs[j],:]
			comp = '{0}-{1}'.format(groupIDs[i], groupIDs[j])
			n1 = t1.shape[0]
			n2 = t2.shape[0]
			deltaI = np.zeros((n1*n2, t1.shape[1]))
			k = 0
			while k < n1*n2:
				for idx1 in range(n1):
					for idx2 in range(n2):
						deltaI[k,:] = brayWrap(t1[idx1,:], t2[idx2,:])
						k += 1
			spMeans =np.round(deltaI.mean(axis=0), 2)
			spSds = np.round(deltaI.std(axis=0, ddof=1), 2)
			spRat = np.round(spMeans/spSds, 2)
			spPct = np.round(spMeans/spMeans.sum()*100, 2)
			tempDF = DataFrame({'sp_mean': spMeans,
						'sp_sd': spSds,
						'ratio': spRat,
						'sp_pct': spPct}, 
						index=[[comp]*len(spNames), spNames])
			tempDF.sort(columns=['sp_pct'], inplace=True, ascending=False)
			tempDF['cumulative'] = np.cumsum(tempDF['sp_pct'])
			tempDF = tempDF[['sp_mean', 'sp_sd', 'ratio', 'sp_pct', 'cumulative']]
			if i==0 and j==i+1:
				finalDF = tempDF
			else:
				finalDF = finalDF.append(tempDF)
			print(comp)
			j += 1
		i += 1
	return finalDF



def brayWrap(x,y):
	temp = np.array([x, y])
	sums = np.sum(temp)
	deltas = np.apply_along_axis(brayFunc, 0, temp, sums=sums)
	return deltas

def brayFunc(m, sums):
	delta = np.abs(m[0]-m[1])/sums
	return 100*delta


	
