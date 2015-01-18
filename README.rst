EcoPy: Python for Ecological Data Analyses
******************************************
**EcoPy** provides tools for ecological data analyses. In general, it focuses on multivariate data analysis, which can be useful in any field, but with particular attention to those methods widely used in ecology. `The homepage, with full documentation and examples, can be found here <http://ecologicalpython.wordpress.com/>`_

What's New
=======
0.0.5
-------
- poca class for princple coordinate analysis
- MDS class for multidimensional scaling (uses isotonic regression from scikit-learn)
- small changes and fixes to previous functions

License
=====
**EcoPy** is distributed under the GNU GPL

Version
=====
0.0.4 - Under development

Examples
======
Transforming a site x species matrix, dividing by site totals::

	import pandas.rpy.common as com
	import ecopy as ep
	varespec = com.load_data('varespec', 'vegan')
	newMat = ep.transform(varespec, method='total', axis=1)

Calculating Bray-Curtis dissimilarities on the new matrix::

	brayMat = ep.distance(newMat, method='bray')

PCA on US Arrests data::
	
	USArrests = com.load_data('USArrests')
	prcomp = ep.pca(USArrests, scaled = True)
	prcomp.summary()
	prcomp.biplot(scale = 0)
	prcomp.biplot(scale = 1, obsNames = True)

Full online documentation is a work in progress

TO-DO
====
- PCoA (MDS)
- RDA
- CCA (INCL. DETRENDED)
- nMDS
- ANOSIM
- SIMPER
- MINIMUM SPANNING TREE (PRIMM's)
- PROCRUSTES ROTATION
- LINEAR/SURFACE ENVIRONMENTAL FITTING
- SPECIES POOLS (ACCUMULATION CURVES)
- MAXENT WRAPPER
- MANY MANY THINGS
