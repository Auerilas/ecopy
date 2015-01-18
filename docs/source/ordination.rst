Ordination
=========

Ecopy contains numerous methods for ordination, that is, plotting points in reduced space. Techniques include, but are not limited to, principle components analysis (PCA), correspondence analysis (CA), principle coordinates analysis (PCoA), and multidimensional scaling (nMDS).

.. py:class:: pca(x, scale=True, varNames=None)

	Takes an input matrix and performs principle components analysis. It will accept either pandas.DataFrames or numpy.ndarrays.  It returns on object of class 'pca', with several methods and attributes. This function uses eigenanalysis of covariance matrices rather than SVD decomposition. NOTE: PCA will NOT work with missing observations, as it is up to the user to decide how best to deal with those.

	**Parameters**

	x: a numpy.ndarray or pandas.DataFrame
		A matrix for ordination, where objects are rows and descriptors/variables as columns. Can be either a pandas.DataFrame or numpy. ndarray

	scale: [True | False]
		Whether or not the columns should be standardized prior to PCA. If 'True', the PCA then operates on a correlation matrix, which is appropriate if variables are on different measurement scales. If variables are on the same scale, use 'False' to have PCA operate on the covariance matrix.

	varNames: list
		If using a numpy.ndarray, pass a list of column names for to help make PCA output easier to interpret. Column names should be in order of the columns in the matrix. Otherwise, column names are represented as integers during summary.

	**Attributes**

	.. py:attribute:: evals
		
		Eigenvalues in order of largest to smallest
		
	.. py:attribute:: evecs
		
		Normalized eigenvectors corresponding to each eigenvalue (i.e. the principle axes)

	.. py:attribute:: scores
		
		Principle component scores of each object (row) on each principle axis. This returns the raw scores :math:`\mathbf{F}` calculated as :math:`\mathbf{F} = \mathbf{YU}` where :math:`\mathbf{U}` is the matrix of eigenvectors and :math:`\mathbf{Y}` are the original observations.

	**Methods**

	.. py:method:: summary_imp()

		Returns a data frame containing information about the principle axes.

	.. py:method:: summary_rot()

		Returns a data frame containing information on axes rotations (i.e. the eigenvectors).

	.. py:method:: summary_corr()

		 Returns a data frame containing the correlation of each variable (column) with each principle axis. For example, the correlation of variable *i* with axis *k* is calculated as :math:`r_{ik} = u_{ik} \sqrt{\lambda_k} / \sqrt{s_i^2}` where :math:`\lambda_k` is the eigenvalue (i.e. variance) associated with axis *k* and :math:`s_i^2` is the variance of variable *i*.

	.. py:method:: summary_desc()

		Returns a data frame containing the cumulative variance explained for each predictor along each principle axis

	.. py:method:: biplot(xax=1, yax=2, type='distance', obsNames=False)

		Create a biplot using a specified transformation.

		xax: integer
			Specifies which PC axis to plot on the x-axis

		yax: integer 
			Specifies which PC axis to plot on the y-axis

		type: ['distance' | 'correlation']
			Type 'distance' plots the raw scores :math:`\mathbf{F}` and the raw vectors :math:`\mathbf{U}` of the first two principle axes. 

			Type 'correlation' plots scores and vectors scaled by the eigenvalues corresponding to each axis: :math:`\mathbf{F\Lambda}^{-0.5}` and :math:`\mathbf{U\Lambda}^{0.5}`, where :math:`\mathbf{\Lambda}` is a diagonal matrix containing the eigenvalues.

		obsNames: [True | False]
			Denotes whether to plot a scatterplot of points (False) or to actually show the names of the observations, as taken from the DataFrame index (True).

	**Examples**

	Principle components analysis of the USArrests data. First, load the data from R using pandas::

		import ecopy as ep
		import pandas.rpy.common as com
		USArrests = com.load_dataset('USArrests')

	Next, run the PCA::

		arrests_PCA = ep.pca(USArrests, scale=True)

	Check the importance of the different axes by examining the standard deviations, which are the square root of the eigenvalues, and the proportions of variance explained by each axis::

		impPC = arrests_PCA.summary_imp()
		print impPC
		            PC1     PC2       PC3     PC4
		Std Dev 1.574878 0.994869 0.597129 0.416449
		Proportion 0.620060 0.247441 0.089141 0.043358
		Cum Prop 0.620060 0.867502 0.956642 1.000000

	Next, examine the eigenvectors and loadings to determine which variables contribute to which axes::

		rotPC = arrests_PCA.summary_rot()
		print rotPC
		         PC1       PC2     PC3        PC4
		Murder 0.535899 0.418181 -0.341233 0.649228
		Assault 0.583184 0.187986 -0.268148 -0.743407
		UrbanPop 0.278191 -0.872806 -0.378016 0.133878
		Rape 0.543432 -0.167319 0.817778 0.089024

	Although the loadings are informative, showing the correlations of each variable with each axis might ease interpretation::

		print arrests_PCA.summary_corr()
		           PC1      PC2      PC3     PC4
		Murder 0.843976 0.658584 -0.537400 1.022455
		Assault 0.580192 0.187021 -0.266773 -0.739593
		UrbanPop 0.166116 -0.521178 -0.225724 0.079942
		Rape 0.226312 -0.069680 0.340563 0.037074

	Then, look to see how much of the variance among predictors is explained by the first two axes::

		print arrests_PCA.summary_desc()
		           PC1      PC2     PC3  PC4
		Murder 0.712296 0.885382 0.926900 1
		Assault 0.843538 0.878515 0.904153 1
		Urban Pop 0.191946 0.945940 0.996892 1
		Rape 0.732461 0.760170 0.998626 1

	Show the biplot using the 'correlation' scaling. Instead of just a scatterplot, use obsNames=True to show the actual names of observations::

		arrests_PCA.biplot(type='correlation', obsNames=True)

	.. figure::  images/corrpca.png
		:align:   center