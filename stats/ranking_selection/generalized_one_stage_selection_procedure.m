% Input: 
% N_vec - vector with sample size
% Output: 
% PCS_D - how many time did we get correct selection
function PCS_vec = generalized_one_stage_selection_procedure(N_vec, k, s, Delta, sigma, iters)

p_stop = 0.9; % PCS high enough 
mu_vec = zeros(k, 1); mu_vec(1:s) = Delta; % create lfc 

num_N = length(N_vec);
PCS_vec = zeros(num_N, 1);
for j=1:N_vec % loop on sample size
    for i=1:iters % loop
        X = randn(k,1) .* (sigma / sqrt(N_vec(j))) + mu_vec;
        [~, I] = sort(X, 'descend');
        if(max(I(1:s)) == s) % correct selection
            PCS_vec(j) = PCS_vec(j)+1;
        end
    end    
    if(PCS_vec(j) > (iters * p_stop)) % here it is high enough, can stop to save time
       stop_at_N = N_vec(j) 
       PCS_vec((j+1):end) = PCS_vec(j); 
       break; 
    end
end
PCS_vec = PCS_vec ./ iters; % normalize



