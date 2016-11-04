# Point-wise Mutual Information (PMI)

This package includes MATLAB implementations for calculating Point-wise Mutual Information (PMI).

Point-wise Mutual Information (PMI)
---------------
**pmi.m** is a program for calculating PMI within a topic or topics. When A is a nonnegative matrix,

    pmi(A, Wtopk_idx, mcnt, epsilon)

returns a matrix of PMI values.
The four parameters are mandatory, and the details are as follows.

* `A`      -  Input matrix (m x n)
* `Wtopk_idx` -  Cell array containing matrix/matrices with indices for top keywords 
* `mcnt`  -  Number of models (e.g. Standard NMF, Sparse NMF, etc.)
* `epsilon`           -  Value for preventing log(0)

One usage example is provided in **example_pmi.m**.

References
----------
1. J. Kim and H. Park. Sparse nonnegative matrix factorization for
  clustering. 2008.
2. D. Newman, J. H. Lau, K. Grieser, and T. Baldwin. Automatic evaluation
   of topic coherence. In Proc. the Annual Conference of the North
   American Chapter of the Association for Computational Linguistics
   (NAACL HLT), pages 100–108, 2010.

Feedback
--------

Please send bug reports, comments, or questions to [Sangho Suh](mailto:sh31659@gmail.com).
Contributions and extentions with new algorithms are welcome.
