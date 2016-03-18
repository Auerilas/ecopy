Version History
==========
0.1.0
^^^^^
- Updated the diversity function
- div_partition function for calculating alpha, beta, and gamma diversity
- spatial_median function for calculation multivariate medians
- fixed a bug in MDS function that provided incorrect results when using monotone transformation
- beta_dispersion function for assessing homogeneity of variances of distance matrices

0.0.9
^^^^^
- Missing data imputation
- nls Python 3 compatibility
- Gower's Euclidean distance for missing data
- ord_plot function for convex hull and line plots of ordination results
- Fully incorporated non-linear regression, including documentation
- Incorporated partial Mantel test in Mantel class
- Global tests of RDA significance
- Updated CCA to include correspondence analysis of residual (unconstrained) variance
- Global tests of CCA significance

0.0.8
^^^^^
- Updated PCA to use SVD instead of eigen decomposition

0.0.7
^^^^^
- CCor
- CCA
- RDA
- RLQ analysis
- Hill and Smith ordination
- weighted mean, variance, scaling


0.0.6
^^^^^
- procrustes test of matrix associations
- anosim class for analysis of similarity
- mantel class for Mantel tests
- corner4 class for fourth corner analysis
- load_data function for importing datasets

0.0.5
^^^^^
- poca class for princple coordinate analysis
- MDS class for multidimensional scaling (uses isotonic regression from scikit-learn)
- small changes and fixes to previous functions

0.0.4
^^^^^
- ca class for simple correspondance analysis

0.0.3
^^^^^
- diversity function for calculation species diversity
- rarefy function for rarefaction

0.0.2
^^^^^
- distance function for calculating distance matrices using a wide variety of coefficients and metrics
- transform function for transforming matrices

0.0.1
^^^^^
- nls class for non-linear regression
- pca class for principle components analysis