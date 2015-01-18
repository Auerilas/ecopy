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
		

