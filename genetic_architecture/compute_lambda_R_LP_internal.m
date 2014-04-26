% Compute sibling relative risk for the LP disease model 
% Function not vectorized
% 
% Input: 
% N - total number of liabilities
% K - number of liabilities which need to be turned 'on'
% h_x - heritability in each liability 
% h_shared_env - shared environment in each liability
% k_R - IBD sharing between two relatives
% mu_l - 'prevalence' of each liability
% mu - prevalence
% 
% Output: 
% lambda_R - risk to relatives
% 
function lambda_R = compute_lambda_R_LP_internal(N, K, h_x, h_shared_env, k_R, mu_l, mu)

if(k_R==0) % singular case - no IBD sharing and no increased risk 
    lambda_R=1;
    return;
end
Z_tab = zeros(2);
Z_tab(2,2) = mu_l^2 * heritability_to_familial_risk(h_x+h_shared_env./k_R, ...
    'liability', mu_l, k_R); % why do we divide here by K_R? since in the function it is multiplied by k_R 
Z_tab(1,2) = mu_l - Z_tab(2,2);  
Z_tab(2,1) = mu_l - Z_tab(2,2);
Z_tab(1,1) = 1 - 2*mu_l + Z_tab(2,2);

lambda_R = 0; % compute joint disease probability for two relatives
for j=K:N % number of 'on' liabilities in person
    for k=K:N % number of 'on' liabilities in relative
        for l=max(0,k+j-N):min(j,k) % number of shared 'on' liabilities
            lambda_R = lambda_R + ...
                Z_tab(2,2)^l * Z_tab(1,2)^(k+j-2*l) * Z_tab(1,1)^(N-k-j+l) * ...
                multinomial(N, [l, k-l, j-l, N-k-j+l]);
        end
    end
end
lambda_R = lambda_R ./ mu^2; % normalize by prevalence squared


