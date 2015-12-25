import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as py

class ca(object):
	'''
	Docstring for function ecopy.ca
	====================
	Conducts correspondance analysis (CA). User supplies 
		an observation x descriptor matrix.

	Use
	----
	ca(x, siteNames=None, spNames=None, scaling=1)

	Returns an object of class ca

	Parameters
	----------
	x:  Data for ordination. Should either be a pandas DataFrame or numpy.ndarray.
		Observations as rows and descriptors as columns. Only
		positive numbers and 0's allowed.
	siteNames: A list of site names
	spNames: A list of species names
	scaling: What type of biplot to produce. See online documentation

	Attributes (see online documentation for descriptions)
	---------
	w_col: column weights of the transformed matrix
	w_row: row weights of the transformed matrix
	evals: eigenvalues of the QQ matrix
	U: column eigenvectors
	Uhat: row eigenvectors
	cumDesc_Sp: The proportion of variance for each species explained by each
		correspondance axis
	cumDesc_Site: The proportion of variance for each site explained by each
		correspondance axis
	siteScores: Site scores along each CA axis
	spScores: Species scores along each CA axis

	Methods
	--------
	summary(): provides a pandas.DataFrame summary table of CA axes
	biplot(coords=False, xax=1, yax=2, type=1, showSp=True, showSite=True, spCol='r', siteCol='k', spSize=12, siteSize=12, xlim=None, ylim=None):
		Produces a biplot of the given CA axes.
		showSp: Whether species should be plotted
		showSite: Whether site should be plotted
		spCol: Color of species text
		siteCol: Color of site text
		spSize: Size of species text
		siteSize: Size of site text
		xlim: Provide a xlim list to override default limits
		ylim: Provide a ylim list to override default limits
		coords: Should the plotting coordinates be returned
		xax: Integer specifying CA Axes to be plotted on the x-axis (Defaults to 1)
		yax: Integer specifying CA Axes to be plotted on the y-axis (Defaults to 2)

	Example
	--------
	import ecopy as ep

	BCI = ep.load_data('BCI')
	bci_ca = ep.ca(BCI)
	print(bci_ca.summary())
	bci_ca.biplot()
	'''
	def __init__(self, x, siteNames=None, spNames=None, scaling=1):
		# if the data is not a dataframe or array, raise error
		if not isinstance(x, (DataFrame, np.ndarray)):
			msg = 'Data must either be pandas.DataFrame or nump.ndarray'
			raise ValueError(msg)
		#  if x is a DataFrame
		if isinstance(x, DataFrame):
			# check NAs
			if x.isnull().any().any():
				msg = 'DataFrame contains null values'
				raise ValueError(msg)
			# check for non-numeric
			if (x.dtypes == 'object').any():
				msg = 'DataFrame can only contain numeric values'
				raise ValueError(msg)
			# convert to a numpy array
			y = np.array(x)
		# if x is array, simple re-assign	
		if isinstance(x, np.ndarray):
			if np.isnan(x).any():
				msg = 'Array contains null values'
				raise ValueError(msg)
			y = x
		# check for negative values
		if y.any() < 0:
			msg ='Matrix cannot contain negative values'
			raise ValueError(msg)
		if scaling not in [1,2]:
			msg = 'type parameter must be 1 or 2'
			raise ValueError(msg)
		if y.shape[0] < y.shape[1]:
			y = y.T
			self.Trans = True
		else:
			self.Trans = False
		pMat = y.astype('float')/y.sum()
		self.w_row = pMat.sum(axis=1)
		self.w_col = pMat.sum(axis=0)
		w_rowA = self.w_row[:,np.newaxis]
		w_colA = self.w_col[np.newaxis,:]
		Q = (pMat - w_rowA*w_colA)/np.sqrt(w_rowA*w_colA)
		self.evals, self.U = np.linalg.eig(Q.T.dot(Q))
		idx = self.evals.argsort()[::-1]
		self.evals = self.evals[idx]
		self.U = self.U[:,idx]
		self.Uhat = Q.dot(self.U).dot(np.diag(self.evals**-0.5))
		self.evals = self.evals[:-1]
		self.U = self.U[:,:-1]
		self.Uhat = self.Uhat[:,:-1]
		if isinstance(x, DataFrame):
			self.siteLabs = x.index
			self.spLabs = x.columns
		else:
			self.siteLabs = ['Site ' + str(x) for x in range(y.shape[0])]
			self.spLabs = ['Sp ' + str(x) for x in range(y.shape[1])]
		if siteNames is not None:
			self.siteLabs = siteNames
		if spNames is not None:
			self.spLabs = spNames
		U2 = self.U.dot(np.diag(self.evals**0.5))
		Uhat2 = self.Uhat.dot(np.diag(self.evals**0.5))
		if self.Trans:
			self.cumDesc_Sp = DataFrame(np.apply_along_axis(lambda x: np.cumsum(x**2) / np.sum(x**2), 1, Uhat2))
			self.cumDesc_Site = DataFrame(np.apply_along_axis(lambda x: np.cumsum(x**2) / np.sum(x**2), 1, U2))
		else:
			self.cumDesc_Sp = DataFrame(np.apply_along_axis(lambda x: np.cumsum(x**2) / np.sum(x**2), 1, U2))
			self.cumDesc_Site = DataFrame(np.apply_along_axis(lambda x: np.cumsum(x**2) / np.sum(x**2), 1, Uhat2))
		if isinstance(x, DataFrame):
			self.cumDesc_Sp.index = x.columns
			self.cumDesc_Site.index = x.index
		self.cumDesc_Sp.columns = ['CA Axis ' + str(x) for x in range(1, len(self.evals) + 1)]
		self.cumDesc_Site.columns = ['CA Axis ' + str(x) for x in range(1, len(self.evals) + 1)]
		V = np.diag(self.w_col**-0.5).dot(self.U)
		Vhat = np.diag(self.w_row**-0.5).dot(self.Uhat)
		F = Vhat.dot(np.diag(self.evals**0.5))
		Fhat = V.dot(np.diag(self.evals**0.5))
		if self.Trans:
			siteCent = Fhat
			spCent = F
			siteOut = V
			spOut = Vhat
			if scaling==1:
				self.siteScores = DataFrame(siteCent, index=self.siteLabs)
				self.spScores = DataFrame(spOut, index=self.spLabs)
			elif scaling==2:
				self.siteScores = DataFrame(siteOut, columns=self.siteLabs)
				self.spScores = DataFrame(spCent, columns=self.spLabs)
		else:
			siteCent = F
			spCent = Fhat
			siteOut = Vhat
			spOut = V
			if scaling==1:
				self.siteScores = DataFrame(siteCent, index=self.siteLabs)
				self.spScores = DataFrame(spOut, index=self.spLabs)
			elif scaling==2:
				self.siteScores = DataFrame(siteOut, index=self.siteLabs)
				self.spScores = DataFrame(spCent, index=self.spLabs)


	def summary(self):
		sds = np.sqrt(self.evals)
		props = self.evals / np.sum(self.evals)
		cumSums = np.cumsum(self.evals) / np.sum(self.evals)
		colNames = ['CA Axis ' + str(x) for x in range(1, len(self.evals)+1)]
		sumTable = DataFrame(np.vstack((sds, props, cumSums)), index=['Inertia', 'Prop.', 'Cum. Prop.'])
		sumTable.columns = colNames
		return sumTable
		
	def biplot(self, xax=1, yax=2, showSp=True, showSite=True, spCol='r', siteCol='k', spSize=12, siteSize=12, xlim=None, ylim=None):
		f, ax = py.subplots()
		if showSite:
			ax.plot(self.siteScores.iloc[:,xax-1], self.siteScores.iloc[:,yax-1], 'ko', ms=0)
			[ax.text(x, y, s, fontsize=siteSize, color=siteCol, ha='center', va='center') for x,y,s in zip(self.siteScores.iloc[:,xax-1], self.siteScores.iloc[:,yax-1], self.siteLabs)]
		if showSp:
			ax.plot(self.spScores.iloc[:,xax-1], self.spScores.iloc[:,yax-1], 'k^', ms=0)
			[ax.text(x,y,s, fontsize=spSize, color=spCol, ha='center', va='center') for x,y,s in zip(self.spScores.iloc[:,xax-1], self.spScores.iloc[:,yax-1], self.spLabs)]
			xmax = max(np.amax(self.siteScores.iloc[:,xax-1]), np.amax(self.spScores.iloc[:,xax-1]))
			xmin = min(np.amin(self.siteScores.iloc[:,xax-1]), np.amin(self.spScores.iloc[:,xax-1]))
			ymax = max(np.amax(self.siteScores.iloc[:,yax-1]), np.amax(self.spScores.iloc[:,yax-1]))
			ymin = min(np.min(self.siteScores.iloc[:,yax-1]), np.min(self.spScores.iloc[:,yax-1]))
		ax.set_xlim([xmin*1.15, xmax*1.15])
		ax.set_ylim([ymin*1.15, ymax*1.15])
		if xlim is not None:
			if not isinstance(xlim, list):
				msg = "xlim must be a list"
				raise ValueError(msg)
			ax.set_xlim(xlim)
		if ylim is not None:
			if not isinstance(ylim, list):
				msg = 'ylim must be a list'
				raise ValueError(msg)
			ax.set_ylim(ylim)
		ax.set_xlabel('CA Axis {!s}'.format(xax))
		ax.set_ylabel('CA Axis {!s}'.format(yax))
		py.show()




	