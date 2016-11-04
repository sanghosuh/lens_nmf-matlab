
clear;
load('example_data.mat')

pmi_vals = zeros(size(Wtopk_idx{1},2),mcnt);

pmi_vals = pmi(A, Wtopk_idx, mcnt, 1e-3);