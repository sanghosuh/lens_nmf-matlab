% Total Document Coverage
% 
% Written by Sangho Suh (sh31659@gmail.com)
%            Dept. of Computer Science and Engineering,
%            Korea University
% 
% Reference: 
% 
% [1] Sangho Suh, Jaegul Choo, Joonseok Lee and Chandan K. Reddy. 
%     L-EnsNMF: Boosted Local Topic Discovery via Ensemble of Nonnegative Matrix Factorization.
% 
% Please send bug reports, comments, or questions to Sangho Suh.
% This comes with no guarantee or warranty of any kind.
%
% Last modified 11/04/2016
%
% this function returns a single value
% i.e. an average value from a number of pmi values
%
% <Inputs>
% 
%       A : [matrix] Input matrix (m x n)
%       term_idx : [matrix] Wtopk matrix with corresponding indices
%       min_nterm : [scalar] min number of keywords
% 
% <Output>
% 
%       qualtopic : [vector] number of documents covered by each topic
%       tot_covrg : [scalar] total number of documents covered by k topics
% 
function [qualtopic,tot_covrg] = compute_total_doc_cvrg(A, term_idx, min_nterm)

        % create 1 x k vector, where k is number of topics
        qualtopic = zeros(1,size(term_idx,2));              
        
        % create k x n matrix, where n is number of documents
        qualtopic_mat = zeros(size(term_idx,2), size(A,2)); 
                             
    % repeat k number of times (where k is number of topics)
    for i=1:size(term_idx,2) 
        
        % get a scalar value that tells how many doc contain min_nterm(c2) number of keywords in a given single topic
        qualtopic(i) = sum( sum( A(term_idx(:,i),:)~=0 )>=min_nterm );  % ... cf.[1]
        
        % get a list of 0's and 1's where 1 indicates that corresponding column has number of keywords more than min_nterm
        qualtopic_mat(i,:) = sum( A(term_idx(:,i),:)~=0 ) >=min_nterm;  % ... cf.[2]        
        
    end
    
    % calculate how many documents are covered by all topics
    tot_covrg = sum( sum(qualtopic_mat) ~=0 ) / size(A,2); % ... cf.[3]
end

        % Breakdown of [1] 
        % (1) term_idx(:,i) is e.g., 3, 32, 47, ... ,132 (where each number refers to index)
        % (2) A(term_idx(:,i),:) is then e.g., A(3,:), A(32,:), A(47,:), ... , (132,:) 
        % (3) A(term_idx(:,i),:)~=0 is then e.g., 
        %                                        [1 0 1 1 0 1 0 0 0 ... ] 
        %                                        [0 1 0 0 1 1 0 0 1 ... ] 
        %                                                ... 
        %     where 1 means it appears at least once in given doc
        %
        % (4) sum(A(term_idx(:,i),:)~=0) is a summation of the above stacked list to give single vector
        %                                  e.g., [4 2 10 3 0 0 2 ...] where each number
        %     
        %     represents the number of keywords each doc contains. So, above
        %     case means doc1 has 4 keywords of given topic, doc2 has 2 keywords of given topic, and so on.
        %
        % (5) sum(sum(A(term_idx(:,i),:)~=0)>=min_nterm) is a summation of the above single vector to give a scalar value,
        %     e.g. 350 where the number represents the total number of
        %     documents that has the same or more number of keywords than min_nterm
        %     thus, qualtopic(1) has the total number of documents covered by Topic 1
        %     So, it is a measure of how good each topic is.
        %     e.g. qualtopic(1) = 24, means, Topic 1's keywords appear in 24 documents in total.
        
        % Breakdown of [2]
        % qualtopic_mat is a matrix where each row represents documents per topic
        %     e.g. qualtopic_mat(1,:) is a 1 x n vector 
        %     (containing 1's and 0's where 1 indicates that corresponding doc(column) has at least min_nterm keywords)
    