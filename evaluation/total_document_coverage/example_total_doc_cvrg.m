
clear;
load('example_data.mat');

min_nterm_list = 1:10; % min number of keywords doc MUST contain (c2 in the paper)

% initialize
qualtopic = {}; totcvrg= {};
qualtopic_mat = zeros(length(min_nterm_list), mcnt);
totcvrg_mat = zeros(length(min_nterm_list), mcnt);
cnt = 0;

for min_nterm = min_nterm_list(:)'
    for idx=1:5
        [qualtopic{idx}, totcvrg{idx}] = compute_total_doc_cvrg(A, Wtopk_idx{idx}, min_nterm);
    end
    
    cnt = cnt + 1;
    
    qualtopic_mat(min_nterm,:) = mean(cell2mat(qualtopic')'); % ... cf.[1] 
    totcvrg_mat(min_nterm,:) = cell2mat(totcvrg);             % ... cf.[2]  
end

% Breakdown of [1]
% e.g. qualtopic{1} is 1 x k where k is number of methods
%      qualtopic{1} = [48 70 2 235 ... 212] where each number indicates total number of documents covered by each method
%      thus, cell2mat(qualtopic') gives a matrix ( k x mcnt) where mcnt is number of methods
%      applying mean() to this matrix gives a vector (1 x mcnt) with each
%      value providing average number of documents covered by (all the topics of) each method 
%      in other words, qualtopic_mat(1,:) gives a vector (1 x k) where each
%      value indicates total number of documents covered by each method
%      when min_nterm is 1
% 
% Breakdown of [2]
%      totcvrg is a vector (1 x mcnt) where mcnt is number of methods
%      e.g. totcvrg{1} is a scalar value and each value of totcvrg
%      represents total number of documents covered by each method
%      
%      
%      
% qualtopic_mat is a calculation of how many documents each topic covered for a given min number of keywords