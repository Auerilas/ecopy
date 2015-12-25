EcoPy: Python for Ecological Data Analyses
******************************************

.. image:: https://zenodo.org/badge/17555/Auerilas/ecopy.svg
   :target: https://zenodo.org/badge/latestdoi/17555/Auerilas/ecopy
   
**EcoPy** provides tools for ecological data analyses. In general, it focuses on multivariate data analysis, which can be useful in any field, but with particular attention to those methods widely used in ecology. `The homepage, with full documentation and examples, can be found here <http://ecopy.readthedocs.org>`_

Install via 'pip install ecopy'

What's New
=======
0.0.9
-------
- missing data imputation methods
- nls Python 3 compatibility
- Gower's Euclidean distance for missing data
- ord_plot function for convex hull or line plots of ordination results
- Fully incorporated non-linear regression, including documentation
- Added the ability to conduct partial Mantel tests in Mantel class
- Global tests of RDA significance
- Updated CCA to include correspondence analysis of residual (unconstrained) variance
- Global tests of CCA significance

License
=====
**EcoPy** is distributed under the MIT license

Version
=====
0.0.9

Examples
======
Transforming a site x species matrix, dividing by site totals::

	import ecopy as ep
	varespec = ep.load_data('varespec')
	newMat = ep.transform(varespec, method='total', axis=1)

Calculating Bray-Curtis dissimilarities on the new matrix::

	brayMat = ep.distance(newMat, method='bray')

PCA on US Arrests data::
	
	USArrests = ep.load_data('USArrests')
	prcomp = ep.pca(USArrests, scaled = True)
	prcomp.summary()
	prcomp.biplot(scale = 0)
	prcomp.biplot(scale = 1, obsNames = True)

Full online documentation is a work in progress

TO-DO
====
- MINIMUM SPANNING TREE
- PROCRUSTES ROTATION
- LINEAR/SURFACE ENVIRONMENTAL FITTING
- MAXENT WRAPPER
- MANY MANY OTHER THINGS
