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
% Last modified 09/26/2016
%
% <Inputs>
%
%        A : Input matrix 
%        k : Number of topics
%        topk : Number of keywords
%        total : Number of iterations
%
% <Outputs>
%
%        Ws, Hs: Results of rank-2 NMF
%        Drs : Set of stage-wise rows with cosine similarity values 
%        Dcs : Set of stage-wise columns with cosine similarity values 
%        As : Updated input matrix
% 
% <Usage Example> 
%
% [Ws, Hs, Drs, Dcs, As] = lens_nmf(A, k, topk, total); 

function [Ws, Hs, Drs, Dcs, As] = lens_nmf(A, k, topk, total)


    % apply l2-normalization and get row-wise and column-wise cosine similarity values 

    A_l2norm_row = bsxfun(@rdivide,A',sqrt(sum((A').^2)))';
    A_cossim_row = A_l2norm_row*A_l2norm_row'; 
      
    A_l2norm_col = bsxfun(@rdivide,A,sqrt(sum(A.^2)));
    A_cossim_col = A_l2norm_col'*A_l2norm_col;
    
    
%%    
    % initialization
    Ws = cell(total, 1); 
    Hs = cell(total, 1); 
    Rs = cell(total, 1);

    As = A; Rs{1} = A;

	vec_norm = 2.0;
	normW = true;
	anls_alg = @anls_entry_rank2_precompute;
	tol = 1e-4;
	maxiter = 10000;

    params_r2 = [];
    params_r2.vec_norm = vec_norm;
    params_r2.normW = normW;
    params_r2.anls_alg = anls_alg;
    params_r2.tol = tol;
    params_r2.maxiter = maxiter;

%%
    for iter=1:(total) % loop for given number of iterations
        
        if iter == 1 
            
            row_idx = datasample(1:size(A,1),1,'Replace',false);
            
        else

            % sample with weight to get one row index
            row_idx = datasample(1:size(A,1),1,'Replace',false,'Weights', sum(cell2mat(Ws'),2));
    
        end            
        
        % get one column index
        col_idx = datasample(1:size(A,2),1,'Replace',false,'Weights', full(sum(abs(Rs{iter}))));

        % update Drs, Dcs with cosine similarity
        [Drs{iter}, Dcs{iter}] = getWeight(Rs{iter},A_cossim_row,A_cossim_col,row_idx,col_idx);

        % update A matrix using Drs,Dcs
        As = update(Rs{iter}, Drs{iter}, Dcs{iter}); 

        [Ws{iter}, Hs{iter}] = nmfsh_comb_rank2(As, rand(size(As,1),2), rand(2,size(As,2)),params_r2);

        if iter <= total
         
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

%%
function [newA] = update(A, Dr, Dc)

    newA = bsxfun(@times,A,Dr);  % multiply each row of A with Dr row
    newA = bsxfun(@times,newA',Dc)';  % multiply each column of A with Dc column  
    
end

%%
function [newA] = update_res_matrix(A, W, H)

    newA = A - W*H;         % get residual matrix   
    newA (newA<0) = 0;      % set any negative element to zero

end

