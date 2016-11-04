% NOTE: This version is NOT in use as it is outdated 
%       (It is before rank-2 NMF has been applied)
% 
% Localizd Ensemble of Nonnegative Matrix Factorization (L-EnsNMF)
%
% Written by Sangho Suh (sh31659@gmail.com)
%            Dept. of Computer Science and Engineering,
%            Korea University
%
% Reference:
%
%  [1] Sangho Suh et al.
%      Boosted L-EnsNMF: Local Topic Discovery via Ensemble of Nonnegative Matrix Factorization.
%      IEEE International Conference on Data Mining 2016.
%
%  [2] Da Kuang Haesun Park
%      Fast Rank-2 Nonnegative Matrix Factorization for Hierarchical Document Clustering
%      International conference on Knowledge Discovery and Data mining 2013
%
% Please send bug reports, comments, or questions to Sangho Suh.
% This code comes with no guarantee or warranty of any kind.
%
% Last modified 11/04/2016
%
% <Inputs>
%
%        A : Input matrix 
%        k : Number of topics
%        topk : Number of keywords
%        iter : Number of iterations
%
% <Outputs>
%
%        Ws, Hs: Results of standard NMF (rank-2 nmf is NOT applied here)
%        Drs : Set of stage-wise rows with cosine similarity values 
%        Dcs : Set of stage-wise columns with cosine similarity values 
%        As : Updated input matrix
% 
% <Usage Example> 
%
% [Ws, Hs, Drs, Dcs, As] = lens_nmf(A, k, topk, iter); 

function [Ws, Hs, Drs, Dcs, As] = lens_nmf(A, k, topk, iter)
   
    params = inputParser;
    params.addParamValue('method','absolute',@(x) ischar(x) );
    par = params.Results;

    % apply l2-normalization and get row-wise and column-wise cosine similarity values 
    
    A_l2norm_row = bsxfun(@rdivide,A',sqrt(sum((A').^2)))';  
    A_cossim_row = A_l2norm_row*A_l2norm_row';               
    
    A_l2norm_col = bsxfun(@rdivide,A,sqrt(sum(A.^2)));       
    A_cossim_col = A_l2norm_col'*A_l2norm_col;               
    
    Ws = cell(iter, 1); 
    Hs = cell(iter, 1); 
    Rs = cell(iter, 1);
    
    As = A; Rs{1} = A;
    

    for iter=1:(iter)
        
        if iter == 1 

            row_idx = datasample(1:size(A,1),1,'Replace',false);
            
        else
            
            % sample with weight to get one row index
            row_idx = datasample(1:size(A,1),1,'Replace',false,'Weights', full(sum(abs(Rs{iter}),2))); 
            
        end                
      
        % get one column index
        col_idx = datasample(1:size(A,2),1,'Replace',false,'Weights', full(sum(abs(Rs{iter}))));   

        % update Drs, Dcs with cosine similarity
        [Drs{iter}, Dcs{iter}] = getWeight(Rs{iter},A_cossim_row,A_cossim_col,row_idx,col_idx);

        % update A matrix using Drs,Dcs
        As = update(Rs{iter}, Drs{iter}, Dcs{iter}); 

        
        if exist('method','var') & strcmp(method,'hals') % if 'method' variable exists and it is 'hals', then do the following
            [Ws{iter}, Hs{iter}] = nmf(As, k,'verbose',0,'method','hals');
        else 
            [Ws{iter}, Hs{iter}] = nmf(As, k,'verbose',0);
        end
        
        if iter <= iter     
            
            % fix W and use unweighted version of A to get H
            [Hs{iter},temp,suc_H,numChol_H,numEq_H] = nnlsm_activeset(Ws{iter}'*Ws{iter},Ws{iter}'*A,0,1,bsxfun(@times,Hs{iter}',1./Dcs{iter})');
            
            % update residual matrix 
            Rs{iter+1} = update_res_matrix(Rs{iter}, Ws{iter},Hs{iter});
            
        end
        
    end
    
end

%% 

function [Dr_new, Dc_new] = getWeight(A,A_cossim_row,A_cossim_col,trm_idx,doc_idx)

    row_smooth = .01;
    col_smooth = .01;

    Dr_new = A_cossim_row(:,trm_idx)*(1-row_smooth)+row_smooth;
    Dc_new = A_cossim_col(:,doc_idx)*(1-col_smooth)+col_smooth;
    
end

function [newA] = update(A, Dr, Dc)

    newA = bsxfun(@times,A,Dr);       % multiply each row of A with Dr row
    newA = bsxfun(@times,newA',Dc)';  % multiply each column of A with Dc column  
    
end

function [newA] = update_res_matrix(A, W, H)

    newA = A - W*H;         % get residual matrix   
    newA (newA<0) = 0;      % set any negative element to zero

end
