class nls(object):
	'''
	Docstring for function ecopy.nls
	=================
	Provides a wrapper for scipy.optimize.leastsq. The function should 
		therefore only return the calculated errors (leastsq automatically 
		minimizes the square of these errors)
	Use
	----
	nls(self, func, p0, xdata, ydata)

	Returns an object of class NLS

	Parameters
	----------
	func: 
	p0: a dictionary of initial parameter estimates
	xdata: an array, list, or series of predictor values
	ydata: an array, list, or series of the response values 

	Attributes
	---------
	cov: covariance matrix of parameters
	inits: initial values
	logLik: log-likelihood of the model
	nparm: number of parameters
	parmEsts: parameter estimates
	parmSE: standard error of parameters
	RMSE: root mean square error
	pvals: p-values of parameters
	tvals: t-values of parameters

	Methods
	--------
	AIC(k=2): returns AIC for the model. Takes argument k (2 by default)
	summary(): printed summary of model output

	Example
	-------
	import ecopy as ep
	import numpy as np

	# import the urchin data
	urchins = ep.load_data('urchins')

	# get the X and Y matrices
	Y = np.array(urchins['Respiration'])*24
	X = np.array(urchins[['UrchinMass', 'Temp']])

	# define the least-squares model per the scipy.optimze methods
	def tempMod(params, X, Y):
		a = params[0]
		b = params[1]
		c = params[2]
		mass = X[:,0]
		temp = X[:,1]
		yHat = a*mass**b*temp**c
		err = Y - yHat
		return(err)

	# initial parameter estimates
	p0 = {'a':1, 'b':1, 'c': 1}

	# model
	tMod = ep.nls(tempMod, p0, X, Y)
	tMod.summary()
	tMod.AIC()
	'''
	def __init__(self, func, p0, xdata, ydata):
		import numpy as np
		from scipy.optimize import leastsq
		import scipy.stats as spst
		# Check the data
		if len(xdata) != len(ydata):
			msg = 'The number of observations does not match the number of rows for the predictors'
			raise ValueError(msg)
		# Check for NA's in the predictor or response
		if np.isnan(xdata).any() == True:
			msg = 'Predicor variable (x) contains missing values'
			raise ValueError(msg)
		if np.isnan(ydata).any() == True:
			msg = 'Response variable (y) contains missing values'
			raise ValueError(msg)
		# Check parameter estimates
		if type(p0) != dict:
			msg = "Initial parameter estimates (p0) must be a dictionary of form p0={'a':1, 'b':2, etc}"
			raise ValueError(msg)
		self.func = func
		self.inits = p0.values()
		self.nobs = len(ydata)
		self.nparm= len(self.inits)
		self.parmNames = p0.keys()
		for i in range(len(self.parmNames)):
			if len(self.parmNames[i]) > 5:
				self.parmNames[i] = self.parmNames[i][0:4]
		# Run the model
		self.mod1 = leastsq(self.func, self.inits, args = (xdata, ydata), full_output=1)
		# Get the parameters
		self.parmEsts = np.round(self.mod1[0], 4)
		# Get the Error variance and standard deviation
		self.RSS = np.sum(self.mod1[2]['fvec']**2)
		self.df = self.nobs - self.nparm
		self.MSE = self.RSS / self.df
		self.RMSE = np.sqrt(self.MSE)
		# Get the covariance matrix
		self.cov = self.MSE * self.mod1[1]
		# Get parameter standard errors
		self.parmSE = np.diag(np.sqrt(self.cov))
		# Calculate the t-values
		self.tvals = self.parmEsts/self.parmSE
		# Get p-values
		self.pvals = (1 - spst.t.cdf(np.abs(self.tvals), self.df))*2
		# Get biased variance (MLE) and calculate log-likehood
		self.s2b = self.RSS / self.nobs
		self.logLik = -self.nobs/2. * np.log(2.*np.pi) - self.nobs/2. * np.log(self.s2b) - 1./(2.*self.s2b) * self.RSS
		del(self.mod1)
		del(self.s2b)       
	
	# Get AIC. Add 1 to the df to account for estimation of standard error
	def AIC(self, k=2):
		print('AIC: {0:5.4f}'.format(-2*self.logLik + k*(self.nparm + 1)))

	# Print the summary
	def summary(self):
		print('')
		print('Non-linear least squares')
		print('Model: ' + self.func.func_name)
		print('Parameters:')
		print('{0:^5} {1:^8} {2:^5} {3:^8} {4:^8}'.format(' ','Estimate', 'Std. Error', 't-value', 'P(>|t|)'))
		for i in range( len(self.parmNames) ):
			print('{0:<5} {1:^8} {2:^8.4f} {3:^8.4} {4:^8.4}'.format(self.parmNames[i], self.parmEsts[i], self.parmSE[i], self.tvals[i], self.pvals[i]))
		print('')
		print('Residual Standard Error: {0:5.4f}'.format(self.RMSE))
		print('Df: {0}'.format(self.df))
