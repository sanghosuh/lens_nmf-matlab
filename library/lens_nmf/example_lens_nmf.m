addpath('../library');

% -------------------------------------------------------------------------
% Usage example for L-EnsNMF
% -------------------------------------------------------------------------

m = 500;
n = 400;

A = rand(m,n);
k = 10;
topk = 10;
iter = 20;

% run L-EnsNMF

[Ws, Hs, Drs, Dcs, As] = lens_nmf(A, k, topk, iter);
