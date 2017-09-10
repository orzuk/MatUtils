% Calculate spreading strategy for current step. 
% We determine how much money to put on each vertex and how much on the global publicity
function [k_vec k_global]=update_strategy(E, X, alpha, k_left, ...
					  old_lambdas, old_lambda_tot)

  N_needed = ceil(length(X)*alpha);
  N_left = N_needed - sum(X); % calc number of non-adopters
				% needed for revolution

				% heuristic: use constant part of money left
				% in each turn
  k_total = 0.3 * k_left;

				% heuristic: use constant part of money left
				% in each turn
%  k_total = k_left / N_left;

				% stupid heuristic: divide money equaly between global and local, and
				% choose a non-adopter randomly.
  k_global = k_total / 2;
  k_vec = zeros(1,length(X));

  nonadopters = find(X==0);
  num_nonadopters = length(nonadopters);

  rand_nonadopt = nonadopters(floor(rand()*num_nonadopters)+1);
  k_vec(rand_nonadopt) = k_total - k_global;

%  [temp, index] = max(old_lambdas .* (X==0));
%  k_vec(index) = k_total / 2;

