function [Y w] = tfidf2( X )
% FUNCTION applies TF-IDF weighting to word count vector matrix.
%
%   [Y w] = tfidf2( X );
%
% INPUT :
%   X        - word count vectors (one column = one document)
%
% OUTPUT :
%   Y        - TF-IDF weighted document-term matrix
%   w        - IDF weights (useful to process other documents)
%

% get inverse document frequencies
w = idf( X );

% TF * IDF
Y = tf( X ) .* repmat( w, 1, size(X,2) );


function Y = tf( X )
% SUBFUNCTION computes word frequencies

Y = X ./ repmat( sum(X,1), size(X,1), 1 );
Y( isnan(Y) ) = 0;


function I = idf(X)
% SUBFUNCTION computes inverse document frequencies

% count DF (document frequency, i.e., the number of words in corpus)
nz = sum( ( X > 0 ), 2 );

% compute idf for each document
I = log( size(X,2) ./ (nz(:) + 1) );