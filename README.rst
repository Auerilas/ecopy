EcoPy: Python for Ecological Data Analyses
******************************************

.. image:: https://zenodo.org/badge/17555/Auerilas/ecopy.svg
   :target: https://zenodo.org/badge/latestdoi/17555/Auerilas/ecopy
   
**EcoPy** provides tools for ecological data analyses. In general, it focuses on multivariate data analysis, which can be useful in any field, but with particular attention to those methods widely used in ecology. `The homepage, with full documentation and examples, can be found here <http://ecopy.readthedocs.org>`_

Install via 'pip install ecopy'

What's New
=======
0.1.0
--------
- Updated the diversity function
- div_partition function for calculating alpha, beta, and gamma diversity
- spatial_median function for calculating multivariate medians
- fixed an error in MDS function wherein nMDS was incorrect
- beta_dispersion function for assessing homogeneity of variances for distance matrices

License
=====
**EcoPy** is distributed under the MIT license

Version
=====
0.1.0

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
