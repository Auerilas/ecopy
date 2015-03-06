Matrix Comparisons
================

Ecopy contains several methods for comparing matrices. Some of these are similar to ordination, while others are more traditional null hypothesis testing.

	- :py:class:`Mantel` (Mantel test)
	- :py:class:`anosim` (Analysis of similarity)
	- :py:function:`simper` (Percentage similarity calculations)
	- :py:class:`procrustes_test` (Procrustes test of matrix correlations)
	- :py:class:`corner4` (Fourth corner analysis)

.. py:class:: Mantel(d1, d2, test='pearson', tail='both', nperm=999)

	Takes two distance matrices for a Mantel test. Returns object of class :py:class:`Mantel`. Calculates the cross-product between lower triangle matrices, using either standardized variables or standardized ranks. The test statistics is the cross-product is divided by 

	.. math::

		\frac{n(n-1)}{2} - 1

	where *n* is the number of objects.

	**Parameters**

	d1:  numpy.ndarray or pandas.DataFrame
		First distance matrix.

	d2: numpy.nadarray or pandas.DataFrame
		Second distance matrix.

	test: ['pearson' | 'spearman']
		'pearson' performs Mantel test on standardized variables.
	
		'spearman' performs Mantel test on standardized ranks.

	tail: ['both' | 'greater' | 'lower']
		'greater' tests the one-tailed hypothesis that correlation is greater than predicted. 
 		
 		'lower' tests hypothsis that correlation is lower than predicted. 
 		
 		'both' is a two-tailed test

 	nperm: int
	 	Number of permutations for the test.

	**Attributes**

	.. py:attribute:: r_obs 

		Observed correlation statistic
	
	.. py:attribute:: pval 

		p-value for the given hypothesis

	.. py:attribute:: tail 

		The tested hypothesis

	.. py:attribute:: test

		 Which of the statistics used, 'pearson' or 'spearman'
	
	.. py:attribute:: perm 

		Number of permutations

	**Methods**

	.. py:classmethod:: summary()

		Prints a summary output table

	**Examples**

	Load the data::

		import ecopy as ep
		v1 = ep.load_data('varespec')
		v2 = ep.load_data('varechem')

	Standardize the chemistry variables and calculate distance matrices::

		v2 = v2.apply(lambda x: (x - x.mean())/x.std(), axis=0)
		dist1 = ep.distance(v1, 'bray')
		dist2 = ep.distance(v2, 'euclidean')

	Conduct the Mantel test::

		mant = ep.Mantel(dist1, dist2)
		print mant.summary()

		Pearson Mantel Test
		Hypothesis = both

		Observed r = 0.305	p = 0.004
		999 permutations

.. py:class:: anosim(dist, factor1, factor2=None, nested=False, nperm=999)

	Conducts analysis of similarity (ANOSIM) on a distance matrix given one or two factors (groups). Returns object of :py:class:`anosim`. Calculates the observed R-statistic as

	.. math::

		R = \frac{r_b - r_w}{\frac{n(n-1)}{4}}

	where :math:`r_w` is the average within-group ranked distances, :math:`r_b` is the average between-group ranked distances, and *n* is the number of objects (rows) in the distance matrix. The factor is then randomly permuted and R recalculated to generate a null distribution.

	**Parameters**

	dist:  numpy.ndarray or pandas.DataFrame
		Square-symmetric distance matrix.

	factor1: numpy.nadarray or pandas.Series or pandas.DataFrame
		First factor.

	factor2: numpy.nadarray or pandas.Series or pandas.DataFrame
		Second factor.

	nested: [True | False]
		Whether factor1 is nested within factor2. If False, then factor1 and factor2 are permuted independently. If Tue, then factor1 is permuted only within groupings of factor2.

 	nperm: int
	 	Number of permutations for the test.

	**Attributes**

	.. py:attribute:: r_perm1

		Permuted R-statistics for factor1

	.. py:attribute:: r_perm2

		Permuted R-statistics for factor1

	.. py:attribute:: R_obs1

		Observed R-statistic for factor1

	.. py:attribute:: R_obs2

		Observed R-statistic for factor2
	
	.. py:attribute:: pval 

		List of p-values for factor1 and factor2
	
	.. py:attribute:: perm 

		Number of permutations

	**Methods**

	.. py:classmethod:: summary()

		Prints a summary output table

	.. py:classmethod:: plot()

		Plots a histogram of R values

	**Examples**

	Load the data::

		import ecopy as ep
		data1 = ep.load_data('dune')
		data2 = com.load_data('dune_env')

	Calculate Bray-Curtis dissimilarity on the 'dune' data, save the 'Management' factor as factor1 and generate factor2::

		duneDist = ep.distance(data1, 'bray')
		group1 = data2['Management']
		group2map = {'SF': 'A', 'BF': 'A', 'HF': 'B', 'NM': 'B'}
		group2 = group1.map(group2map)

	Conduct the ANOSIM::

		t1 = ep.anosim(duneDist, group1, group2, nested=True, nperm=9999)
		print t1.summary()

		ANOSIM: Factor 1
		Observed R = 0.299
		p-value = 0.0217
		9999 permutations

		ANOSIM: Factor 2
		Observed R = 0.25
		p-value = 0.497
		9999 permutations

		t1.plot()

	.. figure::  images/anosimFig.png
		:figwidth: 75 %
		:width: 75 %
		:align:   center

.. py:function:: simper(data, factor, spNames=None)

	Conducts a SIMPER (percentage similarity) analysis for a site x species matrix given a grouping factor. Returns a pandas.DataFrame containing all output for each group comparison. Percent similarity for each species is calculated as the mean Bray-Curtis dissimilarity of each species, given by:

	.. math::

		\Delta_i = \frac{|y_{ik} - y_{ij}|}{\sum_i^n (y_{ik} + y_{ij})}

	The denominator is the total number of individuals in both sites, :math:`y_{ik}` is the number of individuals of species *i* in site *k*, and :math:`y_{ij}` is the number of individuals in site *j*. This is performed for every pairwise combination of sites across two groups and then averaged to yield the mean percentage similarity of the species. This function also calculates the standard deviation of the percentage similarity, the signal to noise ratio (mean / sd) such that a higher ratio indicates more consistent difference, the percentage contribution of each species to the overall difference, and the cumulative percentage difference.

	The output is a multi-indexed DataFrame, with the first index providing the comparison and the second index providing the species. The function lists the index comparison names as it progresses for reference

	**Parameters**

	data:  numpy.ndarray or pandas.DataFrame
		A site x species matrix.

	factor: numpy.nadarray or pandas.Series or pandas.DataFrame or list
		Grouping factor.

	spNames: list
		List of species names. If data is a pandas.DataFrame, then spNames is inferred as the column names. If data is a np.ndarray, then spNames is given integer values unless this argument is provided

	**Examples**

	Conduct SIMPER on the ANOSIM data from above::

		import ecopy as ep
		
		data1 = ep.load_data('dune')
		data2 = com.load_data('dune_env')
		group1 = data2['Management']
		fd = ep.simper(np.array(data1), group1, spNames=data1.columns)

		Comparison indices:
		BF-HF
		BF-NM
		BF-SF
		HF-NM
		HF-SF
		NM-SF

		print fd.ix['BF-NM']

		 sp_mean  sp_sd  ratio  sp_pct  cumulative
		Lolipere     9.07   2.64   3.44   12.43       12.43
		Poatriv      5.47   4.46   1.23    7.50       19.93
		Poaprat      5.25   1.81   2.90    7.19       27.12
		Trifrepe     5.13   2.76   1.86    7.03       34.15
		Bromhord     3.97   2.92   1.36    5.44       39.59
		Bracruta     3.57   2.87   1.24    4.89       44.48
		Eleopalu     3.38   3.57   0.95    4.63       49.11
		Agrostol     3.34   3.47   0.96    4.58       53.69
		Achimill     3.32   2.34   1.42    4.55       58.24
		Scorautu     3.14   2.03   1.55    4.30       62.54
		Anthodor     2.81   3.29   0.85    3.85       66.39
		Planlanc     2.73   2.19   1.25    3.74       70.13
		Salirepe     2.68   2.93   0.91    3.67       73.80
		Bellpere     2.35   1.91   1.23    3.22       77.02
		Hyporadi     2.17   2.45   0.89    2.97       79.99
		Ranuflam     2.03   2.28   0.89    2.78       82.77
		Elymrepe     2.00   2.93   0.68    2.74       85.51
		Callcusp     1.78   2.68   0.66    2.44       87.95
		Juncarti     1.77   2.60   0.68    2.43       90.38
		Vicilath     1.58   1.45   1.09    2.17       92.55
		Sagiproc     1.54   1.86   0.83    2.11       94.66
		Airaprae     1.34   1.97   0.68    1.84       96.50
		Comapalu     1.07   1.57   0.68    1.47       97.97
		Alopgeni     1.00   1.46   0.68    1.37       99.34
		Empenigr     0.48   1.11   0.43    0.66      100.00
		Rumeacet     0.00   0.00    NaN    0.00      100.00
		Cirsarve     0.00   0.00    NaN    0.00      100.00
		Chenalbu     0.00   0.00    NaN    0.00      100.00
		Trifprat     0.00   0.00    NaN    0.00      100.00
		Juncbufo     0.00   0.00    NaN    0.00      100.00

.. py:class:: procrustes_test(mat1, mat2, nperm=999)

	Conducts a procrustes test of matrix associations on two raw object x descriptor matrices. Returns an object of class :py:class:`procrustes_test`. First, both matrices are column-centered. Then, each matrix is divided by the square root of its sum-of-squares. The test statistic :math:`m12^2` is calculated as:

	.. math::

		m12^2 = 1 - (Trace \mathbf{W})^2

	:math:`\mathbf{W}` is the diagonal matrix of eigenvalues for :math:`\mathbf{X'Y}`, which are the two transformed matrices. Then, rows of **X** are randomly permuted and the test statistic recalculated. The *p*-value is the the proportion of random test statistics less than the observed statistic.

	**Parameters**

	mat1:  numpy.ndarray or pandas.DataFrame
		A raw object x descriptor (site x species) matrix.

	factor1: numpy.nadarray or pandas.DataFrame
		A raw object x descriptor (site x descriptor) matrix.

	nperm: int
		Number of permutations in the test

	**Attributes**

	.. py:attribute:: m12_obs

		Observed m12 squared statistic
	
	.. py:attribute:: pval 

		p-value
	
	.. py:attribute:: perm 

		Number of permutations

	**Methods**

	.. py:classmethod:: summary()

		Prints a summary output table

	**Examples**

	Load the data and run the Mantel test::

		import ecopy as ep

		d1 = ep.load_data('varespec')
		d2 = ep.load_data('varechem')
		d = ep.procrustes_test(d1, d2)
		print d.summary()

		m12 squared = 0.744
		p = 0.00701

.. py:class:: corner4(mat1, mat2, nperm=999, model=1, test='both', p_adjustment=None)

	Conducts fourth corner analysis examining associations between species traits and environmental variables. Species traits are given in a species x trait matrix **Q**, species abundances given in a site x species matrix **L**, and environmental traits given in a site x environment matrix **R**. The general concept of fourth corner analysis is to find matrix **D**:

	.. math::

	\begin{bmatrix}
		\mathbf{L} & \mathbf{R} \\
		\mathbf{Q'} & \mathbf{D}
	\end{bmatrix}

	In a simple case, **R** and **Q** contain one environmental variable and one species trait. An expanded correspondance matrix is created following Dray and Legendre (2008). The association between **R** and **Q** is the calculated as follows:

	- If both variables are quantitative, then association is described by Pearson's correlation coefficient *r*
	- If both variables are qualitative, then association is described by :math:`\chi^2` from a contingency table (see Dray and Legendre 2008, Legendre and Legendre 2011)
	- If one variable is quantitative but the other is qualitative, then association is described using the *F*-statistic.

	Significance of the statistics is determined using one of four permutation models (see below). 

	If **R** and **Q** contain more than one variable or trait, then the test iterates through all possible environment-trait combinations. The method automatically determines the appropriates statistics, depending on the data types (float=quantitative or object=qualitative). **NOTE**: As of now, this is quite slow if the number of traits and/or environmental variables is large.

	**Parameters**

	R:  pandas.DataFrame
		A site x variable matrix containing environmental variables for each site. pandas.Series NOT allowed.

	L: numpy.nadarray or pandas.DataFrame
		A site x species matrix of either presence/absence or abundance. Only integer values allowed.

	Q: pandas.DataFrame
		A species x trait matrix containing trait measurements for each species. pandas.Series NOT allowed.

	nperm: int
		Number of permutations in the test

	test: ['both', 'greater', 'lower']
		Which tail of the permutation distribution should be tested against the observed statistic

	p_adjustment: [None, 'bonferroni', 'holm', 'fdr']:
		Which adjustment should be used for multiple comparisons. 'bonferroni' uses Bonferronni correction, 'holm' uses the Bonferroni-Holm correction, and 'fdr' uses the False Discovery Rate correction.

	**Methods**

	.. py:classmethod:: summary()

		Returns a pandas.DataFrame of output

	**Examples**

	Run fourth corner analysis on the aviurba data from R's ade4 package::

		import ecopy as ep

		traits = ep.load_data('avi_traits')
		env = ep.load_data('avi_env')
		sp = ep.load_data('avi_sp')

		fourcorn = ep.corner4(env, sp, traits, nperm=99, p_adjustment='fdr')
		results = fourcorn.summary()

		print results[['Comparison','adjusted p-value']]
		                Comparison  adjusted p-value
		0         farms - feed.hab             1.000
		1       farms - feed.strat             1.000
		2         farms - breeding             1.000
		3        farms - migratory             1.000
		4     small.bui - feed.hab             0.322
		5   small.bui - feed.strat             0.580
		6     small.bui - breeding             1.000
		7    small.bui - migratory             0.909
		8      high.bui - feed.hab             0.111
		...	.......			   ....
		41  veg.cover - feed.strat             1.000
		42    veg.cover - breeding             0.033
		43   veg.cover - migratory             1.000
	