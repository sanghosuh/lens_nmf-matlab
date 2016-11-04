function [W, H, iter, grad] = nmfsh_comb_rank2(A, Winit, Hinit, params)
%%NMFSH_COMB_RANK2 - A fast algorithm for rank-2 nonnegative matrix factorization
% [W, H, iter, grad] = nmfsh_comb_rank2(A, Winit, Hinit, params)
%
% Input parameters
% A: m*n data matrix
% Winit: m*2 matrix for initialization of W
% Hinit: 2*n matrix for initialization of H
% params (optional)
% params.vec_norm (default=2): indicates which norm to use for the normalization of W or H,
%                              e.g. vec_norm=2 means Euclidean norm; vec_norm=0 means no normalization.
% params.normW (default=true): true if normalizing columns of W; false if normalizing rows of H
% params.tol (default=1e-4): tolerance parameter for stopping criterion
% params.maxiter (default=10000): maximum number of iteration times
%
% Output parameters
% W, H: result of rank-2 NMF
% iter: number of ANLS iterations actually used
% grad: relative norm of projected gradient, reflecting the stationarity of the solution
%
% Da Kuang, Haesun Park
% Feb 2013

[m, n] = size(A);
W = Winit;
H = Hinit;
if size(W, 2) ~= 2
	error('Wrong size of W!');
end
if size(H, 1) ~= 2
	error('Wrong size of H!')
end

if ~exist('params', 'var')
	vec_norm = 2.0;
	normW = true;
	tol = 1e-4;
	maxiter = 10000;
else
	if isfield(params, 'vec_norm')
		vec_norm = params.vec_norm;
	else
		vec_norm = 2.0;
	end
	if isfield(params, 'normW')
		normW = params.normW;
	else
		normW = true;
	end
	if isfield(params, 'tol')
		tol = params.tol;
	else
		tol = 1e-4;
	end
	if isfield(params, 'maxiter')
		maxiter = params.maxiter;
	else
		maxiter = 10000;
	end
end

left = H * H';
right = A * H';
for iter = 1 : maxiter
	if rank(left) < 2
		fprintf('The matrix H is singular\n');
		W = zeros(m, 2);
		H = zeros(2, n);
		[U, S, V] = svds(A, 1);
		if sum(U) < 0
			U = -U;
			V = -V;
		end
		W(:, 1) = U;
		H(1, :) = V';
		return;
	end
        W = anls_entry_rank2_precompute(left, right, W);
	norms_W = sqrt(sum(W.^2));
	if min(norms_W) < eps
		error('Error: Some column of W is essentially zero');
	end
	W = bsxfun(@times, W, 1./norms_W);
	left = W' * W;
	right = A' * W;
	if rank(left) < 2
		fprintf('The matrix W is singular\n');
		W = zeros(m, 2);
		H = zeros(2, n);
		[U, S, V] = svds(A, 1);
		if sum(U) < 0
			U = -U;
			V = -V;
		end
		W(:, 1) = U;
		H(1, :) = V';
		return;
	end
        H = anls_entry_rank2_precompute(left, right, H')';
        gradH = left * H - right';
	left = H * H';
	right = A * H';
	gradW = W * left - right;
	if iter == 1
	        initgrad = sqrt(norm(gradW(gradW<=0|W>0))^2 + norm(gradH(gradH<=0|H>0))^2);
		continue;
	else
        	projnorm = sqrt(norm(gradW(gradW<=0|W>0))^2 + norm(gradH(gradH<=0|H>0))^2);
	end
	if projnorm < tol * initgrad
		break;
	end
end
grad = projnorm / initgrad;

if vec_norm ~= 0
	if normW
        	norms = sum(W.^vec_norm) .^ (1/vec_norm);
	        W = bsxfun(@rdivide, W, norms);
        	H = bsxfun(@times, H, norms');
	else    
        	norms = sum(H.^vec_norm, 2) .^ (1/vec_norm);
	        W = bsxfun(@times, W, norms');
        	H = bsxfun(@rdivide, H, norms);
	end
end
