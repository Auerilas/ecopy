Species Diversity
==============

EcoPy contains several methods for estimating species diversity.


.. py:function:: diversity(x, method='shannon', breakNA=True)

	Parameters
	
	**x**: numpy.ndarray or pandas.DataFrame (*required*)
		A site x species matrix, where sites are rows and columns are species.

	**method**: ['shannon' | 'simpson' | 'invSimpson' | 'dominance' | 'spRich' | 'even']
		
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

	**breakNA**: [True | False]
		Whether null values should halt the process. If False, then null values are removed from all calculations.