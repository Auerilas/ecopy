Base Functions
============

EcoPy contains several basic functions:

	- :py:func:`wt_mean`
	- :py:func:`wt_var`
	- :py:func:`wt_scale`

.. py:function:: wt_mean(x, wt=None)
	
	Calculates as weighted mean. Returns a float.

	.. math::

		\mu = \frac{\sum x_iw_i}{\sum w_i}

	**Parameters**
	
	x: numpy.ndarray or list
		A vector of input observations

	wt: numpy.ndarray or list
		A vector of weights. If this vector does not sum to 1, this will be transformed internally by dividing each weight by the sum of weights

	**Example**

	Weighted mean::

		import ecopy as ep
		print ep.wt_mean([1,3,5], [1,2,1])

.. py:function:: wt_var(x, wt=None)
	
	Calculates as weighted variance. Returns a float.

	.. math::

		\sigma^2 = \frac{\sum w_i(x_i - \mu_w)^2}{\sum w_i}

	where :math:`\mu_w` is the weighted mean.

	**Parameters**
	
	x: numpy.ndarray or list
		A vector of input observations

	wt: numpy.ndarray or list
		A vector of weights. If this vector does not sum to 1, this will be transformed internally by dividing each weight by the sum of weights

	**Example**

	Weighted variance::

		import ecopy as ep
		print ep.wt_var([1,3,5], [1,2,1])

.. py:function:: wt_var(x, wt=None)
	
	Returns a vector of scaled, weighted observations.

	.. math::

		z = \frac{x-\mu_w}{\sigma_w}

	where :math:`\mu_w` is the weighted mean and :math:`\sigma_w` is weighted standard deviation (the square root of weighted variance).

	**Parameters**
	
	x: numpy.ndarray or list
		A vector of input observations

	wt: numpy.ndarray or list
		A vector of weights. If this vector does not sum to 1, this will be transformed internally by dividing each weight by the sum of weights

	**Example**

	Weighted variance::

		import ecopy as ep
		print ep.wt_scale([1,3,5], [1,2,1])
