Species Diversity
==============

EcoPy contains several methods for estimating species diversity.

.. py:function:: diversity(x, method='shannon', breakNA=True)
	
	Calculate species diversity for every site in a site x species matrix

	**Parameters**
	
	x: numpy.ndarray or pandas.DataFrame (*required*)
		A site x species matrix, where sites are rows and columns are species.

	method: ['shannon' | 'simpson' | 'invSimpson' | 'dominance' | 'spRich' | 'even']
		*shannon*: Calculates Shannon's H
		
		.. math::
		
			H = -\sum_1^k p_k \log p_k

		where :math:`p_k` is the relative abundance of species *k*

		*simpson*: Calculates Simpson's D

		.. math::

			D = 1 - \sum_1^k p_k^2

		*invSimpson*: Inverse of Simpson's D

		*dominance* Dominance index. :math:`\max p_k`

		*spRich*: Species richness (# of non-zero columns)

		*even*: Evenness of a site. Shannon's H divided by log of species richness.

	breakNA: [True | False]
		Whether null values should halt the process. If False, then null values are removed from all calculations.

	**Example**

	Calculate Shannon diversity of the 'varespec' dataset from R::

		import pandas.rpy.common as com
		import ecopy as ep
		varespec = com.load_data('varespec', 'vegan')
		shannonH = ep.diversity(varespec, 'shannon')

.. py:function:: rarefy(x, method='rarefy', size=None, breakNA=True)
	
	Returns either rarefied species richness or draws a rarefaction curve for each row. Rarefied species richness is calculated based on the smallest sample (default) or allows user-specified sample sizes.

	**Parameters**

	x: numpy.ndarray or pandas.DataFrame (*required*)
		A site x species matrix, where sites are rows and columns are species.

	method: ['rarefy' | 'rarecurve']

		*rarefy*: Returns rarefied species richness.

		.. math::

			S = \sum_1^i 1 - \frac{N}{size}

		where *N* is the total number of individuals in the site, :math:`N_i` is the number of individuals of species *i*, and *size* is the sample size for rarefaction

Next