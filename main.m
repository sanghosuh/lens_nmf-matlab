% Topic Modeling Experiment 
% 
% Methods: StandardNMF, Sparse NMF, Orthogonal NMF, LDA, L-EnsNMF 
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
% [2] H. Kim and H. Park. Sparse non-negative matrix factorizations via
%     alternating non-negativity-constrained least squares for microarray data
%     analysis.
%
% [3] H. Kim and H. Park. Nonnegative matrix factorization based on
%     alternating nonnegativity constrained least squares and active set method.
%     SIAM Journal on Matrix Analysis and Applications, 30(2):713?730,
%     2008.
% 
% [4] D. Kuang and H. Park. Fast rank-2 nonnegative matrix factorization
%     for hierarchical document clustering. In Proc. the ACM SIGKDD
%     International Conference on Knowledge Discovery and Data Mining
%     (KDD), pages 739?747, 2013.
%
% [5] D. M. Blei, A. Y. Ng, and M. I. Jordan. Latent dirichlet allocation.
%     Journal of Machine Learning Research (JMLR), 3:993?1022, 2003.
%
% [6] https://github.com/kimjingu/nonnegfac-matlab
% [7] http://www.cc.gatech.edu/hpark/software/nmf bpas.zip
% [8] http://davian.korea.ac.kr/myfiles/list/Codes/orthonmf.zip
% [9] http://psiexp.ss.uci.edu/research/programs data/toolbox.html
% 
% Please send bug reports, comments, or questions to Sangho Suh.
% This comes with no guarantee or warranty of any kind.
%
% Last modified 11/04/2016
%
% 

% clear work space
clear;

% add path to library folder
addpath('./library');
addpath('./library/nmf');
addpath('./library/lens_nmf');
addpath('./library/ramkis');
addpath('./library/topictoolbox');

% add path to dataset folder
addpath('./data');

% add path to evaluation folder
addpath('./evaluation/topic_coherence')
addpath('./evaluation/total_document_coverage')

%%

% Specify number of experiment(s) to perform
loop = 1;

tic

for numOfLoop=1:loop
    
    loop = loop + 1; % increase count
    
    for choice=1:1  % value of choice belongs to [1,5] where the value indicates dataset
        
        if(choice==1)
            close all;
            clearvars -except loop;
            choice = 1;
            dataname = {};
        elseif(choice==2)
            close all;
            clearvars -except loop;
            choice = 2;
            dataname = {};
        elseif(choice==3)
            close all;
            clearvars -except loop;
            choice = 3;
            dataname = {};
        elseif(choice==4)
            close all;
            clearvars -except loop;
            choice = 4;
            dataname = {};
        elseif(choice==5)
            close all;
            clearvars -except loop;
            choice = 5;
            dataname = {};
        end
        
        %% Decide dataset
        
            if(choice==1)
                load twitter_tlvl10_n2000;
                has_dict = 1;
            elseif(choice==2)
                load enron_tdm_n2000;
                has_dict = 1;
            elseif(choice==3)
                load reuters;
                has_dict = 1;
            elseif(choice==4)
                load vis_paper;
                has_dict = 1;
            elseif(choice==5)
                load 20newsgroups;
                has_dict = 1;
            end
            
          %% initialize
            
            fig_strt_idx = 0;            
            
            k_s = 2;   % number of topics per stage in L-EnsNMF
            topk = 10; % number of top keywords to be listed in order within a topic (denoted as c1 in experiment section)
            total = 5; % number of stages in L-EnsNMF
            
            k_std = k_s*total; % number of total topics
            
            markers = '.ox+*sdv^<>ph'; % markers for graphs
            mcnt = 0;                  % method count
            mname = {};                % method name
            speed = {};                % computing time
          
            A_original = A;            
            
          %% standard NMF (1st method)
            
            % l2 normalization
            A_l2norm = bsxfun(@rdivide,A,sqrt(sum(A.^2)));
            % l1 normalization
            A_l1norm = bsxfun(@rdivide,A,sum(A));
            % tf-idf
            A_idf = tfidf2(A);
            % tf-idf & l2 normalization
            A_l2norm_idf = bsxfun(@rdivide,A,sqrt(sum(A_idf.^2)));
                        
            % use l1-norm/l2-norm weighting/tf-idf
            target_A = A_l1norm;
            % target_A = A_l2norm;
            % target_A = A_idf;
            % target_A = A_l2norm_idf;
            
            mcnt = mcnt + 1; mname{mcnt} = 'StandardNMF'
            
            tic
            [W{mcnt},H{mcnt}] = nmf(target_A, k_std);
            speed{mcnt} = toc;
            
            if ~exist('dict_new','var') % if there is no variable called 'dict_new', then execute the following
                dict_new = strtrim(mat2cell (dictionary, ones (size(dictionary,1),1) , size(dictionary,2) ) );
            end
            
            if exist('titles','var') % if variable called 'titles' exist then execute the following
                % create empty cell
                title_new = cell(length(titles),1);
                % remove spaces within string
                for i=1:length(titles)
                    title_new{i} = strtrim(titles{i}(2,:));
                end
            end
            
          %% sparse NMF (2nd method)
            
            mcnt = mcnt + 1; mname{mcnt} = 'SparseNMF'
            param = [-1 .5];
                        
            tic
            [H{mcnt},W{mcnt}]=nmfsh_comb(target_A', k_std, param); 
            speed{mcnt} = toc;
            
            H{mcnt} = H{mcnt}'; % store transposed version of H
            W{mcnt} = W{mcnt}'; % store transposed version of W

            [Wtopk{mcnt},Htopk{mcnt},DocTopk{mcnt},Wtopk_idx{mcnt}] = parsenmf(W{mcnt},H{mcnt},dict_new,topk);

          %% weak ortho NMF (3rd method)
            
            mcnt = mcnt + 1; mname{mcnt} = 'OrthoNMF' 
            
            Winit=rand(size(A,1),k_std); % create a matrix the size of 'size(A,1) x k_std' as W
            Winit=Winit./repmat(sqrt(sum(Winit.^2,1)),size(A,1),1);  % normalize
            Hinit=rand(k_std,size(A,2)); % create a matrix the size of 'k_std x size(A,1)' as H
            
            beta = 1e-5;
            tic
            [W{mcnt},H{mcnt}]=weakorthonmf(target_A,rand(size(target_A,1),k_std),rand(k_std,size(target_A,2)),k_std,1e-8)
            speed{mcnt} = toc;
            
            [Wtopk{mcnt},Htopk{mcnt},DocTopk{mcnt},Wtopk_idx{mcnt}] = parsenmf(W{mcnt},H{mcnt},dict_new,topk);
            
          %% LDA (4th method)
            
            mcnt = mcnt + 1; mname{mcnt} = 'LDA'
            
            T=k_std;
            % Set the hyperparameters
            BETA=0.01;
            ALPHA=50/T;
            % The number of iterations
            N = 1000;
            % The random seed
            SEED = 3;
            % What output to show (0 = no output; 1 = iterations; 2 = all output)
            OUTPUT = 1;
            
            [ii,jj,ss] = find(A_original);
            
            ss_num = ceil(sum(ss));
            
            ii_new = zeros(ss_num,1);
            jj_new = zeros(ss_num,1);
            
            cnt = 1;
            for i=1:length(ss)
                ii_new(cnt:(cnt+ss(i)-1)) = ii(i);
                jj_new(cnt:(cnt+ss(i)-1)) = jj(i);
                cnt = cnt + ss(i);
            end
            
            tic
            [ W{mcnt},H{mcnt},~ ] = GibbsSamplerLDA( ii_new , jj_new , T , N , ALPHA , BETA , SEED , OUTPUT );
            speed{mcnt} = toc;
            
            H{mcnt} = H{mcnt}';
            
          %% L-EnsNMF (5th method)
            
            mcnt = mcnt + 1;
            
            mname{mcnt} = sprintf('L-EnsNMF')
            
            tic
            [Ws_wgt, Hs_wgt, Drs, Dcs, As] = lens_nmf(target_A, k_s, topk, total);
            speed{mcnt} = toc;
            
            W{mcnt} = []; H{mcnt} = [];
            for i=1:length(Ws_wgt)
                W{mcnt} = [W{mcnt} Ws_wgt{i}];
                H{mcnt} = [H{mcnt}; Hs_wgt{i}];
            end            
            
          %%
            
            Wtopk = {}; Htopk = {}; DocTopk = {}; Wtopk_idx = {};
            topic_num = 1;
            for i=1:mcnt
                [Wtopk{i},Htopk{i},DocTopk{i},Wtopk_idx{i}] = parsenmf(W{i},H{i},dict_new,topk);
                mname{i}
                Wtopk{i}
            end                        
            
          %% Topic Coherence (in PMI)            
            
            % create a zero matrix to store PMI values
            pmi_vals = zeros(size(Wtopk_idx{1},2),mcnt);  
            epsilon = 1e-3    % default value

            for i=1:mcnt
                for topic_idx=1:size(Wtopk_idx{i},2)
                    pmi_vals(topic_idx,i) = compute_pmi_log2(A, Wtopk_idx{i}(:,topic_idx),epsilon);
                end
            end    

          %% Total Document Coverage 
          
            min_nterm_list = 3:10; % min number of keywords doc MUST contain (c2 in the paper)

            qualtopic = {}; totcvrg= {};
            qualtopic_mat = zeros(length(min_nterm_list), mcnt);
            totcvrg_mat = zeros(length(min_nterm_list), mcnt);
            
            % totcvrg_mat is a calculation of how many documents k topics covered
            for min_nterm = min_nterm_list(:)'
                for idx=1:mcnt
                    [qualtopic{idx}, totcvrg{idx}] = compute_total_doc_cvrg(A, Wtopk_idx{idx}, min_nterm);
                end
                qualtopic_mat(min_nterm,:) = mean(cell2mat(qualtopic')');
                totcvrg_mat(min_nterm,:) = cell2mat(totcvrg);
            end

          %% Save Data
%             if(choice==1)
%                 save (sprintf('tlvl10_%03d',loop));
%             elseif(choice==2)
%                 save (sprintf('enron_tdm_n2000_%03d',loop));
%             elseif(choice==3)
%                 save (sprintf('reuters_%03d',loop));
%             elseif(choice==4)
%                 save (sprintf('vis_paper_%03d',loop));
%             elseif(choice==5)
%                 save (sprintf('20newsgroups_%03d',loop));
%             end        
        end        

end    

toc