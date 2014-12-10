import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as py

class ca:
	'''
	Docstring for function ecopy.ca
	====================
	Conducts correspondance analysis (CA). User supplies 
		an observation x descriptor matrix. NOTE: Only
		positive integers and 0 allowed. An error will 
		be raised for non-integer (float) numbers

	Use
	----
	ca(x)

	Returns an object of class ca

	Parameters
	----------
	x:  Data for ordination. Should either be a pandas DataFrame or numpy.ndarray.
		Observations as rows and descriptors as columns. Only
		positive integer data (and zeros) allowed.

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

	Methods
	--------
	summary(): provides a pandas.DataFrame summary table of CA axes
	biplot(coords=False, type=1, showSp=True, showSite=True, spCol='r', siteCol='k', spSize=12, siteSize=12, xlim=None, ylim=None):
		Produces a biplot of the first two CA axes.
		showSp: Whether species should be plotted
		showSite: Whether site should be plotted
		spCol: Color of species text
		siteCol: Color of site text
		spSize: Size of species text
		siteSize: Size of site text
		xlim: Provide a xlim list to override default limits
		ylim: Provide a ylim list to override default limits
		type: What type of biplot to produce. See online documentation
		coords: Should the plotting coordinates be returned

	Example
	--------
	import pandas.rpy.common as com
	import ecopy as ep

	BCI = com.load_data('BCI', 'vegan')
	bci_ca = ep.ca(BCI)
	print bci_ca.summary()
	bci_ca.biplot()
	'''
	def __init__(self, x, scale = True, varNames = None):
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
		# check for integer
		if y.dtype.kind is not 'i':
			msg = 'Matrix contains floats, convert to  integers prior to correspondance analysis'
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
			self.siteLabs = ['Site ' + str(x) for x in range(1, y.shape[0]+1)]
			self.spLabs = ['Sp ' + str(x) for x in range(1, y.shape[1] + 1)]

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

	def summary(self):
		sds = np.sqrt(self.evals)
		props = self.evals / np.sum(self.evals)
		cumSums = np.cumsum(self.evals) / np.sum(self.evals)
		colNames = ['CA Axis ' + str(x) for x in range(1, len(self.evals)+1)]
		sumTable = DataFrame(np.vstack((sds, props, cumSums)), index=['Std. Dev', 'Prop.', 'Cum. Prop.'])
		sumTable.columns = colNames
		return sumTable
		
	def biplot(self, coords=False, type=1, showSp=True, showSite=True, spCol='r', siteCol='k', spSize=12, siteSize=12, xlim=None, ylim=None):
		if type not in [1,2]:
			msg = 'type parameter must be 1 or 2'
			raise ValueError(msg)
		V = np.diag(self.w_col**-0.5).dot(self.U)
		Vhat = np.diag(self.w_row**-0.5).dot(self.Uhat)
		F = Vhat.dot(np.diag(self.evals**0.5))
		Fhat = V.dot(np.diag(self.evals**0.5))
		if self.Trans:
			siteCent = Fhat
			spCent = F
			siteOut = V
			spOut = Vhat
		else:
			siteCent = F
			spCent = Fhat
			siteOut = Vhat
			spOut = V
		f, ax = py.subplots()
		ax.axvline(0, color='k')
		ax.axhline(0, color='k')
		if type==1:
			if showSite:
				ax.plot(siteCent[:,0], siteCent[:,1], 'ko', ms=0)
				[ax.text(x, y, s, fontsize=siteSize, color=siteCol, ha='center', va='center') for x,y,s in zip(siteCent[:,0], siteCent[:,1], self.siteLabs)]
			if showSp:
				ax.plot(spOut[:,0], spOut[:,1], 'k^', ms=0)
				[ax.text(x,y,s, fontsize=spSize, color=spCol, ha='center', va='center') for x,y,s in zip(spOut[:,0], spOut[:,1], self.spLabs)]
				xmax = max(np.amax(siteCent[:,0]), np.amax(spOut[:,0]))
				xmin = min(np.amin(siteCent[:,0]), np.amin(spOut[:,0]))
				ymax = max(np.amax(siteCent[:,1]), np.amax(spOut[:,1]))
				ymin = min(np.min(siteCent[:,1]), np.min(spOut[:,1]))
			else:
				xmin = np.min(siteCent[:,0])
				xmax = np.max(siteCent[:,0])
				ymin = np.min(siteCent[:,1])
				ymax = np.max(siteCent[:,1])
			ax.set_xlim([xmin+0.15*xmin, xmax+0.15*xmax])
			ax.set_ylim([ymin+0.15*ymin, ymax+0.15*ymax])
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
		if type==2:
			if showSp:
				ax.plot(spCent[:,0], spCent[:,1], 'ko', ms=0)
				[ax.text(x, y, s, fontsize=spSize, color=spCol, ha='center', va='center') for x,y,s in zip(spCent[:,0], spCent[:,1], self.spLabs)]
			if showSite:
				ax.plot(siteOut[:,0], siteOut[:,1], 'k^', ms=0)
				[ax.text(x, y, s, fontsize=siteSize, color=siteCol, ha='center', va='center') for x,y,s in zip(siteOut[:,0], siteOut[:,1], self.siteLabs)]
				xmax = max(np.amax(spCent[:,0]), np.amax(siteOut[:,0]))
				xmin = min(np.amin(spCent[:,0]), np.amin(siteOut[:,0]))
				ymax = max(np.amax(spCent[:,1]), np.amax(siteOut[:,1]))
				ymin = min(np.min(spCent[:,1]), np.min(siteOut[:,1]))
			else:
				xmax = np.amax(spCent[:,0])
				xmin = np.amin(spCent[:,0])
				ymax = np.amax(spCent[:,1])
				ymin = np.min(spCent[:,1])
			ax.set_xlim([xmin+0.15*xmin, xmax+0.15*xmax])
			ax.set_ylim([ymin+0.15*ymin, ymax+0.15*ymax])
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
		ax.set_xlabel('CA Axis 1')
		ax.set_ylabel('CA Axis 2')
		py.show()
		if coords:
			return {'F': siteCent, 'Fhat': spCent, 'V': spOut, 'Vhat': siteOut}




	