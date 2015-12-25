Regression
===================

EcoPy provides a wrapper for scipy.optimize.leastsq. This wrapper allows users to specify a non-linear function and receive parameter estimates, statistics, log-likelihoods, and AIC.

	- :py:class:`nls`

.. py:class:: nls(func, p0, xdata, ydata)

	nls takes a function (func), initial parameter estimates (p0), predictor variables (xdata), and a response (ydata) and passes these to scipy.optimize.leastsq. It returns an object of class :py:class:`nls`.

	**Parameters**

	func: function
		A function that returns the quantity to be minimized. See documentation for scipy.optimize.leastsq.

	p0: dictionary
		A dictionary of initial parameter estimates for every parameter to be estimated in the function.

	xdata: numpy.ndarray
		A numpy.ndarray of predictor variables. See example below for how to include multiple predictors.

	ydata: numpy.ndarray
		A numpy.ndarray of the response variable.

	**Attributes**

	.. py:attribute:: cov
		
		Variance-covariance matrix of parameters
		
	.. py:attribute:: inits

		Initial parameter estimates

	.. py:attribute:: logLik

		Log-likelihood of the function

	.. py:attribute:: nparm

		Number of parameters estimated

	.. py:attribute:: parmEsts

		Parameter estimates

	.. py:attribute:: parmSE

		Standard error of each parameter

	.. py:attribute:: RMSE

		Root mean square error

	.. py:attribute:: pvals

		p-values of each parameter

	.. py:attribute:: tvals

		t-values of each parameter

	**Methods**

	.. py:classmethod:: AIC(k=2)

		Returns AIC for the given model. Argument k determines the correction applied to the number of parameters

	.. py:classmethod:: summary

		Returns a regression summary table

	**Examples**

	First, load the urchin data::

		import ecopy as ep
		import numpy as np

		urchins = ep.load_data('urchins')

	Next, make the X and Y matrices::

		Y = np.array(urchins['Respiration])*24
		X = np.array(urchins[['UrchinMass', 'Temp']])

	Define the least-squares function to be optimized::

		def tempMod(params, X, Y):
			a = params[0]
			b = params[1]
			c = params[2]
			mass = X[:,0]
			temp = X[:,1]
			yHat = a*mass**b*temp**c
			err = Y - yHat
			return(err)

	Create a dictionary of initial estimates for each parameter::

		p0 = {'a':1, 'b':1, 'c': 1}

	Run the model and check the summary tables::

		tMod = ep.nls(tempMod, p0, X, Y)
		tMod.summary()

		Non-linear least squares
		Model: tempMod
		Parameters:
		 	Estimate 	Std. Error	t-value			P(>|t|)
		a    	 0.0002	 0.0002	 0.8037	 0.4302
		c    	 0.3346	 0.1485	 2.2533	 0.0345
		b    	 1.5209	 0.3448	 4.4112	 0.0002

		Residual Standard Error:  0.0371
		Df: 22

		tMod.AIC()

		AIC:  -88.9664797962