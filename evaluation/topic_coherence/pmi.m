% Point-wise Mutual Information (PMI)
% 
% Written by Sangho Suh (sh31659@gmail.com)
%            Dept. of Computer Science and Engineering,
%            Korea University
% 
% Reference: 
% 
% [1] J. Kim and H. Park. Sparse nonnegative matrix factorization for
%     clustering. 2008.
% [2] D. Newman, J. H. Lau, K. Grieser, and T. Baldwin. Automatic evaluation
%     of topic coherence. In Proc. the Annual Conference of the North
%     American Chapter of the Association for Computational Linguistics
%     (NAACL HLT), pages 100?108, 2010.
% 
% Please send bug reports, comments, or questions to Sangho Suh.
% This comes with no guarantee or warranty of any kind.
%
% Last modified 10/17/2016
%
% this function returns a single value
% i.e. an average value from a number of pmi values
%
% <Inputs>
% 
%       A : Input matrix (m x n)
%       Wtopk_idx : Cell array containing matrix with indices for top keywords 
%       mcnt : Number of methods
%       epsilon : (Default: 1e-3)    
% 
% <Output>
% 
%       pmi_vals : A matrix of PMI values
%
function [pmi_vals] = pmi(A, Wtopk_idx, mcnt, epsilon)

    % create a zero matrix to store PMI values
    pmi_vals = zeros(size(Wtopk_idx{1},2),mcnt);  

    for i=1:mcnt
        for topic_idx=1:size(Wtopk_idx{i},2)
            pmi_vals(topic_idx,i) = compute_pmi_log2(A, Wtopk_idx{i}(:,topic_idx),epsilon);
        end
    end    
    
end    