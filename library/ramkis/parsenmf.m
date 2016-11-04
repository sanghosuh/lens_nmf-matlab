function [ Wtopk,Htopk,DocTopk,Wtopk_idx ] = parsenmf( W,H,vocab,topk)
%W is of size mxk
%H is of size kxn
k = size(W,2);      % k <- num of topics
Wtopk=cell(topk,k); % create empty matrix of size, topk x k, for storing 
Wtopk_idx = zeros(topk,k);
n = size(H,2);      % n <- num of documents
Htopk=zeros(min([topk size(W,2)]),n); % zeroes of k x 512

if ~iscell(vocab)
    vocab = strtrim(mat2cell (vocab, ones (size(vocab,1),1) , size(vocab,2) ) );
end

%generate topk words based on W
for i=1:k
    [~,idx]=sort(W(:,i),'descend');
    for j=1:topk
%         try
            Wtopk{j,i}=vocab{idx(j)};
%         catch
%             Wtopk{j,i}=vocab(idx(j),:); % for vis_paper, dict_new of which is not cell type but matrix 
%         end
        Wtopk_idx(j,i) = idx(j);
    end
end
%find top k topics for every document/tweet based on H
for j=1:n
    [~,idx]=sort(H(:,j),'descend');
    Htopk(:,j)=idx(1:min([topk size(W,2)]));
end

DocTopk = zeros(topk,k);

for j=1:k
    [~,idx]=sort(H(j,:),'descend');
    DocTopk(:,j) = idx(1:topk);
end

end
