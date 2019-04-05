EcoPy: Python for Ecological Data Analyses
******************************************

.. image:: https://zenodo.org/badge/17555/Auerilas/ecopy.svg
   :target: https://zenodo.org/badge/latestdoi/17555/Auerilas/ecopy
   
**EcoPy** provides tools for ecological data analyses. In general, it focuses on multivariate data analysis, which can be useful in any field, but with particular attention to those methods widely used in ecology. `The homepage, with full documentation and examples, can be found here <http://ecopy.readthedocs.io>`_

Install via 'pip install ecopy'

What's New
=======
I'm back! I apologize for letting EcoPy languish for a few years. I was busy as a post-doc trying to secure my faculty position. Now I've done that, I can get back to running EcoPy and expanding it. I can't promise I will be devoted to the project at all times, but I am now able to spend more time on it. My hope is to recruit help so that EcoPy expands rapidly in the future.

.. 0.1.2.4
.. --------
.. - Recompiled the isotonic regression using updated Cython for compatability with Python 3.7

0.1.2.3
--------
- Fixed compatibility problems in functions cca(), simper(), and transform()
- Added Ochiai's binary coefficient to distance() function
- Added bioenv() function

0.1.2.2
--------
- More Python 3.x compatibility
- Fix typos in code and examples on readthedocs. Thorough code check


License
=====
**EcoPy** is distributed under the MIT license

Version
=====
0.1.2.4

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
	prcomp = ep.pca(USArrests, scale = True)
	prcomp.biplot(type = 'distance')
	prcomp.biplot(type = 'correlation')

Full online documentation is a work in progress

TO-DO
====
Incorporate DECORANA and TWINSPAN into EcoPy
---------------------------------------------

1. I have modified write_cep to handle integer row names (as is common in pandas dataframes)
2. The pre-processing code all works for DECORANA
3. **Need to get decorana fortran function working on UNIX systems (.exe binary only works for Windows)**
4. **Need to get TWINSPAN functional**

Procrusted Rotation
-------------------

Linear/surface environmental fitting
-------------------------------------

MaxEnt Wrapper
--------------

Clustering
----------
