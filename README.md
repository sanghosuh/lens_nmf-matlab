# Localized Ensemble of Nonnegative Matrix Factorization (L-EnsNMF) 

This package includes MATLAB implementations of localized ensemble of Nonnegative Matrix Factorization (L-EnsNMF).

Localized Ensemble of Nonnegative Matrix Factorization (L-EnsNMF)
---------------
**lens_nmf.m** is a program for executing L-EnsNMF algorithm. When A is a nonnegative matrix,

    lens_nmf(A, k, topk, iter)

returns the L-EnsNMF of A with k as a target lower-rank. 
The four parameters are mandatory, and the details are as follows.

* `A`      -  Input matrix (m x n)
* `k` -  Number of topics
* `topk`  -  Number of keywords
* `iter`           -  Number of iterations

One usage example is provided in **example_lens_nmf.m**.

References
----------
1. Sangho Suh, Jaegul Choo, Joonseok Lee and Chandan K. Reddy.
   Boosted L-EnsNMF: Local Topic Discovery via Ensemble of Nonnegative Matrix Factorization.

2. Da Kuang Haesun Park.
   Fast Rank-2 Nonnegative Matrix Factorization for Hierarchical Document Clustering.

3. Jingu Kim, Yunlong He, and Haesun Park.
   Algorithms for Nonnegative Matrix and Tensor Factorizations: A Unified View 
   Based on Block Coordinate Descent Framework.
   Journal of Global Optimization, 58(2), pp. 285-319, 2014.

Feedback
--------

Please send bug reports, comments, or questions to [Sangho Suh](mailto:sh31659@gmail.com).
Contributions and extentions with new algorithms are welcome.
