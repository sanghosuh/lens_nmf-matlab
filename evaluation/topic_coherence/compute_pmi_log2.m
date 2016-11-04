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
%       term_idx_mat : Wtopk matrix with corresponding indices 
%       epsilon : Value to prevent log(0) (Default: 1e-3)    
% 
% <Output>
% 
%       pmi_val : PMI value
%
function [pmi_val] = compute_pmi_log2(A, term_idx_mat,epsilon)

    % create a term-by-document matrix consisting of only 1 and 0
    % where 1 denotes term occurrence in the corresponding column
    A(A~=0)=1;

    % topic coherence for a topic (Keyword-wise similarity)
    if size(term_idx_mat,2)==1 % if a single column, i.e., one topic
        
        term_idx = term_idx_mat(:);
        
        % create an array for storing PMI values between pairs of keywords within a topic
        pmi_vals = zeros( length(term_idx)*(length(term_idx)-1)/2, 1); 
        cnt = 1;
        
        % calculate PMI values on pairs of keywords in a topic
        for i=1:(length(term_idx)-1)
            for j=(i+1):length(term_idx)
                % if keywords do NOT co-occur in document(s)
                if( ( mean ( A(term_idx(i),:) & A(term_idx(j),:) ) ) == 0 ) 
                    pmi_vals(cnt) = 0;
                % if keywords co-occur in document(s), calculate their PMI value
                else 
                    pmi_vals(cnt) = log2( (mean(A(term_idx(i),:) & A(term_idx(j),:))+epsilon) / (mean(A(term_idx(i),:))*mean(A(term_idx(j),:)))   );
                end
                cnt = cnt + 1;
            end
            % get average from array of PMI values
            pmi_val = mean(pmi_vals); 
        end
    
    % topic coherence for topics (Topic-wise similarity)
    else  % if more than one topic, find the similarity between topics
        
        % create an array for storing PMI values between topics
        pmi_val_tlvl = zeros( size(term_idx_mat,2)*(size(term_idx_mat,2)-1)/2, 1 ); 
        cnt_tlvl = 1;
        
        % calculate PMI values between keywords of topics
        for k=1:(size(term_idx_mat,2)-1)
            for l=(k+1):size(term_idx_mat,2)
                term_idx1 = term_idx_mat(:,k);
                term_idx2 = term_idx_mat(:,l);

                % create an array for storing PMI values between keywords of topics
                pmi_vals = zeros( length(term_idx1)*length(term_idx2), 1 );
                cnt = 1;
                for i=1:length(term_idx1)
                    for j=1:length(term_idx2)
                       % if keywords do NOT co-occur in document(s)
                       if( (mean(A(term_idx1(i),:) & A(term_idx2(j),:))) == 0 ) 
                           pmi_vals(cnt) = 0;
                       % if keywords co-occur in document(s), calculate their PMI value
                       else 
                           pmi_vals(cnt) = log2((mean(A(term_idx1(i),:) & A(term_idx2(j),:))+epsilon) / (mean(A(term_idx1(i),:))*mean(A(term_idx2(j),:))));
                       end
                        cnt = cnt + 1;
                    end
                end
                pmi_val_tlvl(cnt_tlvl) = mean(pmi_vals); 
                cnt_tlvl = cnt_tlvl + 1;
            end
        end
        % get average from array of PMI values        
        pmi_val = mean(pmi_val_tlvl);
    end
end