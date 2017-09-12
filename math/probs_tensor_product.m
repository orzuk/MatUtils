% Generate a joint probability distribution from a set of two distributions 
% 
% Input: 
% x_vec - first vector
% x_probs - probabilities for first vector 
% x_vec - first vector
% x_probs - probabilities for first vector 
% T_marginal - table of marginal probabilities 
% 
% Output: 
% z_vec - joint vector 
% z_probs - probabilities of joint vector 
% 
function [z_vec z_probs] = probs_tensor_product(x_vec, x_probs, y_vec, y_probs, T_marginal)

[x_len x_width] = size(x_vec); [y_len y_width] = size(y_vec);
z_vec = zeros(x_len*y_len, x_width+y_width); 
z_vec(:,x_width+1:end) = repmat(y_vec, x_len, 1); 
z_vec(:,1:x_width) = reshape(repmat(x_vec, 1, y_len)', x_width, x_len*y_len)'; 

independent_flag = 1; 
if(exist('T_marginal', 'var') && (~isempty(T_marginal)))
	independent_flag = 0;
end
if(independent_flag) % assume independence
	z_probs = reshape(repmat(x_probs, 1, y_len)', x_len*y_len, 1) .* repmat(y_probs, x_len, 1) ; 
else % dependence via marginals. Here x and y must be of the same length 
	z_probs = ones(x_len*y_len,1);
    num_vals = size(T_marginal, 1);
	for i=1:x_width
		for j=0:num_vals-1
			for k=0:num_vals-1
				cur_inds = find((z_vec(:,i) == j) & (z_vec(:,x_width+i) == k));
				z_probs(cur_inds) = z_probs(cur_inds) .* T_marginal(j+1,k+1,i); 
			end 
		end				
	end	

end
