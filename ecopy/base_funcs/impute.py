import numpy as np
from pandas import DataFrame, Series

def impute(Y, method='mice', m=5, delta=0.0001, niter=100):
	"""
	Docstring for function ecopy.impute
	====================
	Performs univariate missing data imputation using one of several methods described below

	Use
	----
	impute(Y, method='mice', m=5, delta=0.0001, niter=100)

	Returns a numpy.ndarray

	Parameters
	----------
	Y: A pandas.DataFrame, pandas.Series, or numpy.ndarray containing missing data.
	method: Imputation method
		'mean': Replaces all missing values with the mean of the corresponding column
		'median': Replaces all missing values with the median of the corresponding column
		'multi_norm': Calculates the mean and covariance matrix of the fully observed data.
				Replaces missing values with random draws from this multivariate
				normal distribution.
		'univariate': Bayesian imputation of missing observations based on posterior draws
				from parameter estimates using all other variables as predictors. Operates
				independently for for each column. Assumes that X_miss has no NA's, however
				if X_miss has NA's, NA's are replaced with mean values.
		'monotone': Conducts monotone imputation for longitudinally structured missing values
				using Bayesian regression (i.e. method univariate).
		'mice': Conducts the MICE algorithm for data imputation using the 'univariate' method.
	m: Integer denoting number of multiple imputations to perform.
	delta: Ridge parameter to avoid singular matrices
	niter: Number of iterations for the MICE algorithm

	Example
	--------	
	import ecopy as ep
	import numpy as np
	import pandas as pd

	data = ep.load_data('urchins')
	massNA = np.random.randint(0, 24, 5)
	respNA = np.random.randint(0, 24, 7)

	data.loc[massNA, 'UrchinMass'] = np.nan
	data.loc[respNA, 'Respiration'] = np.nan

	# MICE algorithm
	imputedData = ep.impute(data, 'mice') 
	# convert back to pandas dataframes
	imputedFrame = [pd.DataFrame(x, columns=data.columns) for x in imputedData]

	# Replace with means
	meanImpute = ep.impute(data, 'mean')
	"""

	if not isinstance(Y, (Series, DataFrame, np.ndarray)):
		msg = 'Y must be a pandas.Series, pandas.DataFrame, or numpy.ndarray'
		raise ValueError(msg)
	if delta < 0 or delta >= 0.1:
		msg = 'k should be 0 <= k < 0.1'
		raise ValueError(msg)
	if method not in ['mean', 'median', 'multi_norm', 'univariate', 'monotone', 'mice']:
		msg = 'method is not an approved imputation technique'
		raise ValueError(msg)
	fullData = np.array(Y)
	completeObs = np.apply_along_axis(lambda x: np.sum(np.isnan(x)), 1, fullData)
	obsData = fullData[completeObs==0,:]
	if obsData.shape[0]==1:
		msg = 'Only one row of complete observations'
		raise ValueError(msg)
	if method=='mean':
		result = meanFunc(fullData, obsData)
	if method=='median':
		result = medianFunc(fullData, obsData)
	if method=='multi_norm':
		result = []
		k = 0
		while k < m:
			tempR = multinormFunc(fullData, obsData)
			result.append(tempR)
			k += 1
	if method=='univariate':
		result = []
		k = 0
		while k < m:
			tempR = unipostFunc(fullData, delta)
			result.append(tempR)
			k += 1
	if method=='monotone':
		result = []
		k = 0
		while k < m:
			tempR = monotoneFunc(fullData, delta)
			result.append(tempR)
			k += 1
	if method=='mice':
		result = []
		k = 0
		while k < m:
			tempR = miceFunc(fullData, delta, niter)
			result.append(tempR)
			k += 1
	return result


def meanFunc(fullMat, obsMat):
	z = fullMat.copy()
	means = obsMat.mean(axis=0)
	for i in range(len(means)):
		idx = np.isnan(z[:,i])
		z[idx,i] = means[i]
	return z

def medianFunc(fullMat, obsMat):
	z = fullMat.copy()
	medians = np.median(obsMat, axis=0)
	for i in range(len(medians)):
		idx = np.isnan(z[:,i])
		z[idx,i] = medians[i]
	return z

def multinormFunc(fullMat, obsMat):
	z = fullMat.copy()
	means = obsMat.mean(axis=0)
	cov = np.cov(obsMat, rowvar=0)
	randoms = np.random.multivariate_normal(means, cov, fullMat.shape[0]*10)
	for i in range(len(means)):
		idx = np.isnan(z[:,i])
		randVals = np.random.choice(randoms[:,i], idx.sum(), replace=False)
		z[idx,i] = randVals
	return z


def unipostFunc(fullMat, delta):
	z = fullMat.copy()
	r = fullMat.copy()
	for i in range(z.shape[1]):
		yj = z[:,i]
		xj = np.delete(z, i, axis=1)
		idx = np.isnan(yj)
		if idx.sum()==0:
			r[:,i] = yj
		else:
			yj_obs = yj[~idx]
			xj_obs = xj[~idx,:]
			xj_nan = np.apply_along_axis(lambda x: np.sum(np.isnan(x)), 1, xj_obs)
			xj_obs = xj_obs[xj_nan==0,:]
			yj_obs = yj_obs[xj_nan==0]
			S = xj_obs.T.dot(xj_obs)
			V = np.linalg.pinv(S + np.diag(np.diag(S))*delta)
			b_obs = V.dot(xj_obs.T).dot(yj_obs)
			g = np.random.chisquare(xj_obs.shape[0] - xj_obs.shape[1])
			resid = yj_obs - xj_obs.dot(b_obs)
			MSE = np.sqrt((resid.T.dot(resid))/g)
			w1 = np.random.randn(xj_obs.shape[1])
			G = np.linalg.cholesky(V)
			b_star = b_obs + MSE*w1.dot(V.T)
			w2 = np.random.randn(idx.sum())
			xmis = xj[idx]
			for k in range(xmis.shape[1]):
				idx2 = np.isnan(xmis[:,k])
				xmis[idx2,k] = np.nanmean(xmis[:,k])
			y_star = xmis.dot(b_star) + w2*MSE
			yj[idx] = y_star
			r[:,i] = yj
	return r


def monotoneFunc(fullMat, delta):
	z = fullMat.copy()
	r = fullMat.copy()
	nMiss = np.apply_along_axis(lambda x: np.sum(np.isnan(x)), 0, z)
	sortID = nMiss.argsort()
	origID = np.arange(z.shape[1])
	orderedZ = z[:,sortID]
	orderedR = r[:,sortID]
	originalOrder = origID[sortID]
	for i in range(orderedZ.shape[1]):
		yj = orderedZ[:,i]
		idx = np.isnan(yj)
		if i==0:
			yj[idx] = yj.mean()
			orderedR[:,i] = yj
		else:
			xj = orderedR[:,range(i)]
			yj_obs = yj[~idx]
			xj_obs = xj[~idx,:]
			S = xj_obs.T.dot(xj_obs)
			V = np.linalg.pinv(S + np.diag(np.diag(S))*delta)
			b_obs = V.dot(xj_obs.T).dot(yj_obs)
			g = np.random.chisquare(xj_obs.shape[0] - xj_obs.shape[1])
			resid = yj_obs - xj_obs.dot(b_obs)
			MSE = np.sqrt((resid.T.dot(resid))/g)
			w1 = np.random.randn(xj_obs.shape[1])
			G = np.linalg.cholesky(V)
			b_star = b_obs + MSE*w1.dot(V.T)
			w2 = np.random.randn(idx.sum())
			xmis = xj[idx]
			for k in range(xmis.shape[1]):
				idx2 = np.isnan(xmis[:,k])
				xmis[idx2,k] = np.nanmean(xmis[:,k])
			y_star = xmis.dot(b_star) + w2*MSE
			yj[idx] = y_star
			orderedR[:,i] = yj
	rFinal = orderedR[:,originalOrder]
	return rFinal


def miceFunc(fullMat, delta, niter):
	z = fullMat.copy()
	r = fullMat.copy()
	for i in range(z.shape[1]):
		yj = z[:,i]
		idx = np.isnan(yj)
		if idx.sum()==0:
			r[:,i] = yj
		else:
			yobs = yj[~idx]
			yj[idx] = np.random.choice(yobs, idx.sum(), replace=False)
			r[:,i] = yj
	t = 0
	while t < niter:
		for i in range(r.shape[1]):
			yj = r[:,i]
			xj = np.delete(r, i, axis=1)
			idx = np.isnan(yj)
			if idx.sum()==0:
				r[:,i] = yj
			else:
				yj_obs = yj[~idx]
				xj_obs = xj[~idx,:]
				xj_nan = np.apply_along_axis(lambda x: np.sum(np.isnan(x)), 1, xj_obs)
				xj_obs = xj_obs[xj_nan==0,:]
				yj_obs = yj_obs[xj_nan==0]
				S = xj_obs.T.dot(xj_obs)
				V = np.linalg.pinv(S + np.diag(np.diag(S))*delta)
				b_obs = V.dot(xj_obs.T).dot(yj_obs)
				g = np.random.chisquare(xj_obs.shape[0] - xj_obs.shape[1])
				resid = yj_obs - xj_obs.dot(b_obs)
				MSE = np.sqrt((resid.T.dot(resid))/g)
				w1 = np.random.randn(xj_obs.shape[1])
				G = np.linalg.cholesky(V)
				b_star = b_obs + MSE*w1.dot(V.T)
				w2 = np.random.randn(idx.sum())
				xmis = xj[idx]
				for k in range(xmis.shape[1]):
					idx2 = np.isnan(xmis[:,k])
					xmis[idx2,k] = np.nanmean(xmis[:,k])
				y_star = xmis.dot(b_star) + w2*MSE
				yj[idx] = y_star
				r[:,i] = yj
		t += 1
	return r



