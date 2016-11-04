# Total Document Coverage

This package includes MATLAB implementations for calculating total document coverage.

Total Document Coverage
---------------
**total_doc_cvrg.m** is a program for calculating total document coverage of a topic or topics. When A is a nonnegative matrix,

    total_doc_cvrg(A, Wtopk_idx, mcnt, epsilon)

returns two variables:
* 'qualtopic' - 
* 'qualtopic_mat' - 


The four parameters are mandatory, and the details are as follows.

* `A`      -  Input matrix (m x n)
* `Wtopk_idx` -  Cell array containing matrix/matrices with indices for top keywords 
* `mcnt`  -  Number of models (e.g. Standard NMF, Sparse NMF, etc.)
* `epsilon`           -  Value for preventing log(0)

One usage example is provided in **example_pmi.m**.

References
----------
1. Sangho Suh, Jaegul Choo, Joonseok Lee and Chandan K. Reddy. L-EnsNMF: Boosted Local Topic Discovery via Ensemble of Nonnegative Matrix Factorization.

Feedback
--------
Please send bug reports, comments, or questions to [Sangho Suh](mailto:sh31659@gmail.com).
Contributions and extentions with new algorithms are welcome.
