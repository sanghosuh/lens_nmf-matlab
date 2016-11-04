# Localized Ensemble of Nonnegative Matrix Factorization (L-EnsNMF) 
============

This package includes MATLAB implementations for L-EnsNMF as well as other state-of-the-art topic modeling methods.
The methods are as follows. 

* 'Standard NMF'
* 'Sparse NMF' 
* 'Orthogonal NMF'
* 'LDA'
* 'L-EnsNMF'

## Experiment
---------------
**main.m** is a program for running experiment(s) on dataset(s) using the above methods.
The program returns topic keywords of each topic modeling method and evaluation results.
Some of them are as follows.

* 'Wtopk' - cell array (1 x mcnt) where mcnt is the number of methods
* 'speed' - cell array (1 x mcnt) where mcnt is the number of methods
* 'totcvrg_mat' - matrix (k x mcnt) where c1 is the number of keywords

## Evaluation
---------------
This package includes two evaluation measures:
 
* 'Topic coherence'
* 'Total document coverage'

While one can see how it can be used from **main.m**, they have separate,simple usage example files, **example_pmi.m** and **example_total_doc_cvrg.m**.
More details about them can be found [here](https://github.com/sanghosuh/lens_nmf-matlab/evaluation/).

## References
----------
1.   Sangho Suh, Jaegul Choo, Joonseok Lee and Chandan K. Reddy. 
     L-EnsNMF: Boosted Local Topic Discovery via Ensemble of Nonnegative Matrix Factorization.
     International Conference on Data Mining(ICDM), 2016.

2.   D. Kuang and H. Park. Fast rank-2 nonnegative matrix factorization
     for hierarchical document clustering. In Proc. the ACM SIGKDD
     International Conference on Knowledge Discovery and Data Mining
     (KDD), pages 739?747, 2013.

3.   Jingu Kim, Yunlong He, and Haesun Park. Algorithms for Nonnegative Matrix and Tensor Factorizations: 
     A Unified View Based on Block Coordinate Descent Framework. Journal of Global Optimization, 58(2), pp. 285-319, 2014.

4.   J. Kim and H. Park. Sparse nonnegative matrix factorization for
     clustering. 2008.

5.   D. Newman, J. H. Lau, K. Grieser, and T. Baldwin. Automatic evaluation
     of topic coherence. In Proc. the Annual Conference of the North
     American Chapter of the Association for Computational Linguistics
     (NAACL HLT), pages 100–108, 2010.

6.   H. Kim and H. Park. Sparse non-negative matrix factorizations via
     alternating non-negativity-constrained least squares for microarray data
     analysis.

7.   H. Kim and H. Park. Nonnegative matrix factorization based on
     alternating nonnegativity constrained least squares and active set method.
     SIAM Journal on Matrix Analysis and Applications, 30(2):713?730,
     2008.

8.   D. M. Blei, A. Y. Ng, and M. I. Jordan. Latent dirichlet allocation.
     Journal of Machine Learning Research (JMLR), 3:993?1022, 2003.

9.   https://github.com/kimjingu/nonnegfac-matlab
10.  http://www.cc.gatech.edu/~hpark/software/nmf bpas.zip
11.  http://davian.korea.ac.kr/myfiles/list/Codes/orthonmf.zip
12.  http://psiexp.ss.uci.edu/research/programs data/toolbox.html


Feedback
--------
Please send bug reports, comments, or questions to [Sangho Suh](mailto:sh31659@gmail.com).
Contributions and extentions with new algorithms are welcome.
