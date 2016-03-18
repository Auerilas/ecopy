import numpy as np
from pandas import DataFrame 

def distance(x, method='euclidean', transform="1", breakNA=True):
	'''
	Docstring for function ecopy.distance
	========================
	Computes a dissimilarity matrix from a given matrix,
		where distances are calculated among rows. Multiple
		similarity/dissimilarity metrics are possible, and it is
		up to the user to choose the correct one. In the case of similarity
		metrics, dissimilarity is calculated as 1-D unless otherwise specified.
		See Legendre and Legendre for details on coefficients.

	For binary coefficients, A is # shared similarities, B and C are #'s of
		differences, and D is # of shared absences

	Most quantitative (i.e. non-binary) metrics will work with binary data.
		In fact, some methods are redundant. For example, bray = sorensen 
		when bray is passed a binary data frame. However, both methods are 
		provided here for simplicity and clarity

	Use
	----
	distance(x, method='euclidean', transform="1", breakNA=True)

	Returns a numpy.ndarray square, symmetric matrix

	Parameters
	----------
	x:  numpy array or pandas dataframe with observations as rows
		and descriptors as columns
	method: a method used for calculating similarities
		euclidean: calculates euclidean distance between rows (default)
				sqrt(d1^2 + d2^2 + ... + dn^2)
		gow_euclidean: calculates euclidean distance between rows, removing NAs
				sqrt(delta*(x1-x2)**2/sum(delta))
				where delta = 1 if both observations present, 0 otherwise
		chord: calculates the chord distance between rows
				this is the euclidean distance of normalized vectors
		manhattan: manhattan distance between rows
				sum(abs(xik - xjk)) over all k
		meanChar: Czekanowski's mean character difference
				1/M*(sum(abs(xik - xjk))) over all k
				where M is the number of columns
		whittaker: Whittaker's index of association
				0.5*sum(abs(x1-x2)) where x1 and x2 are
				first standardized by vector totals (x1/sum(x1))
		canberra: Canberra metric
				sum(abs(xik-xjk) / (xik+xjk))*1/M 
				where M is the number of species present at both sites
		hellinger: Hellinger distance
				same as the chord but square-root transformed first
		mod_gower: modified Gower distance
				same as meanChar except M is the number of 
				columns that are not double zero
 		simple: simple matching of binary data (SIMILARITY)
				(A + D)/(A + B + C + D)
		rogers: Rogers and Tanimoto coefficient for binary data (SIMILARITY)
				(A + D) / (A + 2B + 2C + D)
		sokal: Sokal and Sneath coefficient for binary data (SIMILARITY)
				(2A + 2D)/(2A + B + C + 2D)
		jaccard: Jaccard's coefficient for binary data (SIMILARITY)
				A/(A + B + C)
		sorensen: Sorensen's coefficient for binary data (SIMILARITY)
				2A/(2A + B + C)
		kulczynski: Kulczynski's coefficient (SIMILARITY)
				0.5*(sum(min(xi,xj))/sum(xi) + sum(min(xi,xj))/sum(xj))
		bray: Bray-Curtis coefficient (SIMILARITY)
				2*sum(min(xi, xj))/(sum(xi) + sum(xj))
		gower: Gower assymetrical coefficient (SIMILARITY)
				(1-(abs(xik - xjk))/(max(xk)-min(xk)))*1/M
				where k is species and M is the number of columns.
				The denominator is the maximum of species k minus the min of 
				species k in the entire matrix. Note that double zero columns
				are excluded here.
		NOTE: All metrics denoted as SIMILARITY are transformed to dissimilarity using
			D = 1 - Similarity
	transform: how should dissimilarities be transformed.
		"1": D (default)
		"sqrt": sqrt(D)
	breakNA: should the process halt if the matrix contains any NAs?
		if False, then NA's undergo pairwise deletion during distance calculation,
		such that when calculating the distance between two rows, if any
		species is missing from a row, that species is removed from both rows

	Example
	--------
	import ecopy as ep
	varespec = ep.load_data('varespec')
	distance(varespec, method='euclidean')

	# for binary data
	varespec[varespec>0] = 1
	distance(varespec, method='jaccard)
	'''
	listofmethods =['euclidean', 'gow_euclidean', 'simple', 'rogers', 'sokal', 'jaccard', 'sorensen', 'kulczynski', 'bray', 'gower', 'chord', 'manhattan', 'meanChar', 'whittaker', 'canberra', 'hellinger', 'mod_gower']
	if not isinstance(breakNA, bool):
		msg = 'removaNA argument must be boolean'
		raise ValueError(msg)
	if method not in listofmethods:
		msg = 'method argument {0!s} is not an accepted metric'.format(method)
		raise ValueError(msg)
	if not isinstance(x, (DataFrame, np.ndarray)):
		msg = 'x argument must be a numpy array or pandas dataframe'
		raise ValueError(msg)
	if isinstance(x, DataFrame):
		if (x.dtypes == 'object').any():
			msg = 'DataFrame can only contain numeric values'
			raise ValueError(msg)
		x = np.array(x)
	if breakNA:
		if np.isnan(np.sum(x)):
			msg = 'Matrix contains NA values'
			raise ValueError(msg)
	if np.min(np.sum(x, axis=1))==0:
		msg = 'One row is entirely zeros, distance calculations will be meaningless'
		raise ValueError(msg)
	if transform not in ['1', 'sqrt']:
		msg = 'transform argument must be "1" or "sqrt"'
		raise ValueError(msg)
	if method in ['simple', 'rogers', 'sokal', 'jaccard', 'sorensen']:
		if np.any((x != 0) & (x != 1)):
			msg = 'For method {0}, data must be binary'.format(method)
			raise ValueError(msg)
	x = x.astype('float')
	if method == 'euclidean':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = eucFunc(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'gow_euclidean':
		distMat = np.zeros((x.shape[0]. x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,:]
				x2 = x[j,:]
				delta = ~(np.isnan(x1) + np.isnan(x2))
				delta = delta.astype(int)
				x1[np.isnan(x1)] = -999
				x2[np.isnan(x2)] = -999
				distMat[i,j] = eucGow(x1, x2, delta, transform)
				distMat[j,i] = distMat[i,j]
	if method == 'simple':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				A, B, C, D = matchMat(x1, x2)
				distMat[i,j] = simpleSim(A, B, C, D, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'rogers':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				A, B, C, D = matchMat(x1, x2)
				distMat[i,j] = rogerSim(A, B, C, D, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'sokal':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				A, B, C, D = matchMat(x1, x2)
				distMat[i,j] = sokalSim(A, B, C, D, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'jaccard':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				A, B, C, D = matchMat(x1, x2)
				distMat[i,j] = jaccardSim(A, B, C, D, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'sorensen':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				A, B, C, D = matchMat(x1, x2)
				distMat[i,j] = sorenSim(A, B, C, D, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'kulczynski':
		if (x<0).any():
			msg = 'Distances are meaningless for negative numbers'
			raise ValueError(msg)
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = kulSim(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'bray':
		if (x<0).any():
			msg = 'Distances are meaningless for negative numbers'
			raise ValueError(msg)
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = braySim(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'gower':
		if (x<0).any():
			msg = 'Distances are meaningless for negative numbers'
			raise ValueError(msg)
		distMat = np.zeros((x.shape[0], x.shape[0]))
		R = np.apply_along_axis(lambda z: np.max(z) - np.min(z), 0, x)
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				R = R[~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				R = R[~np.isnan(x2)]
				distMat[i,j] = gowerSim(x1, x2, R, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'chord':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = chordDis(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'manhattan':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = manDist(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'meanChar':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = charDist(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'whittaker':
		if (x<0).any():
			msg = 'Distances are meaningless for negative numbers'
			raise ValueError(msg)
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = whitDist(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'canberra':
		if (x<0).any():
			msg = 'Distances are meaningless for negative numbers'
			raise ValueError(msg)
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = canDist(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'hellinger':
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = chordDis(np.sqrt(x1), np.sqrt(x2), transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)
	if method == 'mod_gower':
		if (x<0).any():
			msg = 'Distances are meaningless for negative numbers'
			raise ValueError(msg)
		distMat = np.zeros((x.shape[0], x.shape[0]))
		for i in range(0, distMat.shape[0]):
			distMat[i,i] = 0
			for j in range(i+1, distMat.shape[0]):
				x1 = x[i,~np.isnan(x[i,:])]
				x2 = x[j,~np.isnan(x[i,:])]
				x1 = x1[~np.isnan(x2)]
				x2 = x2[~np.isnan(x2)]
				distMat[i,j] = m_gowDist(x1, x2, transform)
				distMat[j,i] = distMat[i,j]
		return(distMat)

def eucFunc(d1, d2, t):
	d = d1-d2
	eucD = np.sqrt(np.sum(np.square(d)))
	if t == "1":
		return eucD
	if t == "sqrt":
		return np.sqrt(eucD)

def eucGow(d1, d2, Delta, t):
	d = d1-d2
	eucGow = np.sqrt(np.sum(Delta*d**2)/Delta.sum())
	if t == "1":
		return eucGow
	if t == "sqrt":
		return np.sqrt(eucGow)

def matchMat(d1, d2):
	A = float(np.sum((d1 == 1) & (d2 == 1)))
	B = float(np.sum((d1 == 1) & (d2 == 0)))
	C = float(np.sum((d1 == 0) & (d2 == 1)))
	D = float(np.sum((d1 == 0) & (d2 == 0)))
	return A, B, C, D

def simpleSim(A, B, C, D, t):
	S = (A + D) / (A + B + C + D)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)

def rogerSim(A, B, C ,D, t):
	S = (A + D) / (A + 2*B + 2*C + D)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)

def sokalSim(A, B, C, D, t):
	S = (2*A + 2*D) / (2*A + B + C + 2*D)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)

def jaccardSim(A, B, C, D, t):
	S = A/(A + B + C)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)

def sorenSim(A, B, C, D, t):
	S = 2*A/(2*A + B + C)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)

def kulSim(d1, d2, t):
	A = np.sum(d1)
	B = np.sum(d2)
	W = np.sum(np.minimum(d1, d2))
	S = 0.5*(W/A + W/B)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)
		
def braySim(d1, d2, t):
	A = np.sum(d1)
	B = np.sum(d2)
	W = np.sum(np.minimum(d1, d2))
	S = (2*W)/(A + B)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)

def gowerSim(d1, d2, R, t):
	diffs = 1 - (np.abs(d1-d2)/R)
	isabs = d1+d2 ==0
	S = np.sum(~isabs*diffs)/np.sum(~isabs)
	if t == "1":
		return 1 - S
	if t == "sqrt":
		return np.sqrt(1 - S)

def chordDis(d1, d2, t):
	norm1 = 1/(np.sqrt(np.sum(d1**2))) * d1
	norm2 = 1/(np.sqrt(np.sum(d2**2))) * d2
	d = norm1-norm2
	chordD = np.sqrt(np.sum(np.square(d)))
	if t == "1":
		return chordD
	if t == "sqrt":
		return np.sqrt(chordD)

def manDist(d1, d2, t):
	manD = np.sum(np.abs(d1 - d2))
	if t == "1":
		return manD
	if t == "sqrt":
		return np.sqrt(manD)

def charDist(d1, d2, t):
	charD = 1./len(d1) * np.sum(np.abs(d1 - d2))
	if t == "1":
		return charD
	if t == "sqrt":
		return np.sqrt(charD)

def whitDist(d1, d2, t):
	d1 = d1/np.sum(d1)
	d2 = d2/np.sum(d2)
	whitD = 0.5*sum(np.abs(d1-d2))
	if t == "1":
		return whitD
	if t == "sqrt":
		return np.sqrt(whitD)

def canDist(d1, d2, t):
	isabs = d1+d2==0
	d1 = d1[~isabs]
	d2 = d2[~isabs]
	canD = np.sum(np.abs(d1-d2)/(d1+d2)) * 1./np.sum(~isabs)
	if t == "1":
		return canD
	if t == "sqrt":
		return np.sqrt(canD)

def m_gowDist(d1, d2, t):
	isabs = d1+d2==0
	charD = 1./np.sum(~isabs) * np.sum(np.abs(d1 - d2))
	if t == "1":
		return charD
	if t == "sqrt":
		return np.sqrt(charD)