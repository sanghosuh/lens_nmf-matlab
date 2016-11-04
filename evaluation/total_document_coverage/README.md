# Total Document Coverage

This package includes MATLAB implementations for calculating total document coverage.

Total Document Coverage
----------------
**compute_total_doc_cvrg.m** is a program for calculating total document coverage of a topic or topics.

    compute_total_doc_cvrg(A, term_idx, min_nterm)

returns two values, which are as follows.

* `qualtopic` - Number of documents covered by each topic (1 x k) where k is number of topics
* `tot_cvrg` - Total number of documents covered by k topics (1 x 1)

The three parameters are mandatory, and the details are as follows.

* `A`      -  Input matrix (m x n)
* `term_idx` -  Cell array containing matrix/matrices with indices for top keywords 
* `min_nterm`  -  Minimum number of keywords document must have (c2 in References[1])

One usage example is provided in **example_total_doc_cvrg.m**.

References
----------
1. Sangho Suh, Jaegul Choo, Joonseok Lee and Chandan K. Reddy. L-EnsNMF: Boosted Local Topic Discovery via Ensemble of Nonnegative Matrix Factorization. International Conference on Data Mining (ICDM), 2016.

Feedback
--------
Please send bug reports, comments, or questions to [Sangho Suh](mailto:sh31659@gmail.com).
Contributions and extentions with new algorithms are welcome.
