% Compute joint probs from marginals
% 
% Input: 
% x_vec - set of variables vectors 
% f_vec - vector of marginal probabilities
% normalize_flag - normalize output probs. to sum to one (default: NO)  
%
% Output: 
% p_x_vec - probability of each value assuming Pr(x) = \prod_i [f(i)^(x(i)) + (1-f(i))*(1-x(i))] 
%
function p_x_vec = x_to_prob_vec(x_vec, f_vec, normalize_flag)

num_x = size(x_vec, 1);
p_x_vec = repmat(f_vec, num_x, 1); % probability of each of 2^N possibilities (according to independence model)
p_x_vec = prod(p_x_vec .^ x_vec,2) .* prod((1-p_x_vec) .^ (1-x_vec),2);

if(~exist('normalize_flag', 'var') || isempty(normalize_flag))
    normalize_flag = 0; 
end
if(normalize_flag)
    p_x_vec = p_x_vec ./ sum(p_x_vec);
end
