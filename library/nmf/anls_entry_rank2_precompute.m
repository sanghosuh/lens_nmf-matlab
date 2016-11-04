function H = anls_entry_rank2_precompute(left, right, H)

% left: 2*2
% right: n*2
% Returning H of size n*2 also
%
% Da Kuang, Haesun Park
% Feb 2013

n = size(right, 1);

solve_either = zeros(n, 2);
solve_either(:, 1) = right(:, 1) * (1./left(1,1));
solve_either(:, 2) = right(:, 2) * (1./left(2,2));
cosine_either = bsxfun(@times, solve_either, [sqrt(left(1,1)), sqrt(left(2,2))]);
choose_first = (cosine_either(:, 1) >= cosine_either(:, 2));
solve_either(choose_first, 2) = 0;
solve_either(~choose_first, 1) = 0;

%H = (left \ right')';
if abs(left(1,1)) < eps & abs(left(1,2)) < eps
	error('Error: The 2x2 matrix is close to singular or the input data matrix has tiny values');
else
	if abs(left(1,1)) >= abs(left(1,2)) % pivot == 1, extracting cosine
		t = left(2,1) / left(1,1);
		a2 = left(1,1) + t * left(2,1);
		b2 = left(1,2) + t * left(2,2);
		d2 = left(2,2) - t * left(1,2);
		if abs(d2/a2) < eps % a2 is guaranteed to be positive
			error('Error: The 2x2 matrix is close to singular');
		end
		e2 = right(:, 1) + t * right(:, 2);
		f2 = right(:, 2) - t * right(:, 1);
	else % pivot == 2, extracting sine
		ct = left(1,1) / left(2,1);
		a2 = left(2,1) + ct * left(1,1);
		b2 = left(2,2) + ct * left(1,2);
		d2 = -left(1,2) + ct * left(2,2);
		if abs(d2/a2) < eps % a2 is guaranteed to be positive
			error('Error: The 2x2 matrix is close to singular');
		end
		e2 = right(:, 2) + ct * right(:, 1);
		f2 = -right(:, 1) + ct * right(:, 2);
	end
	H(:, 2) = f2 * (1/d2);
	H(:, 1) = (e2 - b2 * H(:, 2)) * (1/a2);
end

use_either = ~all(H>0, 2);
H(use_either, :) = solve_either(use_either, :);
