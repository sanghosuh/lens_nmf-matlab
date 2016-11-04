# Evaluation Measures for Topic Modeling Algorithms

This package includes MATLAB implementations for two evaluation measures:

* `Topic coherence (PMI)`
* `Total document coverage`

Topic Coherence
----------------
**pmi.m** is a program for calculating PMI within a topic or topics. When A is a nonnegative matrix,

    pmi(A, Wtopk_idx, mcnt, epsilon)

returns a matrix of PMI values.
The four parameters are mandatory, and the details are as follows.

* `A`      -  Input matrix (m x n)
* `Wtopk_idx` -  Cell array containing matrix/matrices with indices for top keywords 
* `mcnt`  -  Number of models (e.g. Standard NMF, Sparse NMF, etc.)
* `epsilon`           -  Value for preventing log(0)

One usage example is provided in **example_pmi.m**.


Total Document Coverage
----------------
**compute_total_doc_cvrg.m** is a program for calculating total document coverage for a topic or topics.

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
