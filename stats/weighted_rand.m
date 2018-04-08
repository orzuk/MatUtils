% Produce random numbers drawn according to some weights. 
% We assume weights are positive and indices are drawn from [1:L] where L is the number of
% different weights (weights can be a matrix)
%
% Input:
% weights - probability weights. L is length 
% n - number of random values to produce
%
% Output:
% r - vector of random integers from 1:L according to weights 
%
function r = weighted_rand(weights, n)

if(issparse(weights)) % a faster version for sparse weights vectors 
    [~, J, S] = find(weights);
    r = J(weighted_rand(S, n)); 
    return;
end
    
L = size(weights, 2); m = size(weights, 1);
cum_weights = cumsum(weights,2);
cum_weights = cum_weights ./ repmat(cum_weights(:,end), 1, L); % normalize weights

r = zeros(m, n); % initilize indices
tmp_r = rand(m, n); % randomize
if(L == 1) % only one weight (binary search doesn't work) ..
    r(:) = 1;
else
    for w_ind = 1:m
        if(cum_weights(w_ind,1)==1) % take first, Binary search doesn't work 
            r(w_ind,:)=1;
        else
            r(w_ind,:) = bsearch(cum_weights(w_ind,:), tmp_r(w_ind,:));
            f = find(tmp_r(w_ind,:) > cum_weights(w_ind,r(w_ind,:)));
            r(w_ind, f) = r(w_ind, f) + 1;
        end
    end
end

% Old slow way:
% for iter = 1:n
%     tmp_ind = bsearch(cum_weights, tmp_r)
%     for i=L:-1:1
%         r(tmp_r(:,iter) < cum_weights(:,i), iter) = i;
%     end
% end

