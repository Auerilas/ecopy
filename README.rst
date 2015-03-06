EcoPy: Python for Ecological Data Analyses
******************************************
**EcoPy** provides tools for ecological data analyses. In general, it focuses on multivariate data analysis, which can be useful in any field, but with particular attention to those methods widely used in ecology. `The homepage, with full documentation and examples, can be found here <http://ecopy.readthedocs.org>`_

What's New
=======
0.0.6
-------
- procrustes_test for procrustes test of matrix associations
- load_data function for loading data
- anosim class for analysis of similarity
- mantel class for Mantel tests

License
=====
**EcoPy** is distributed under the GNU GPL

Version
=====
0.0.6 - Under development

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
- RDA
- CCA (INCL. DETRENDED)
- MINIMUM SPANNING TREE (PRIMM's)
- PROCRUSTES ROTATION
- LINEAR/SURFACE ENVIRONMENTAL FITTING
- MAXENT WRAPPER
- MANY MANY OTHER THINGS
