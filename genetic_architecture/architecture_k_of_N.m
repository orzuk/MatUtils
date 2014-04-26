% Compute all statistics for a special architecture which requires K of N
% architectures to be 'on'. Computation is analytic.
%
% Input:
% N - number of loci
% K - minimal number of 'on' loci needed
% mu_l - 'prevalence' of each input 'disease'
%
% Output:
% lambda_R - familial relative risk for the combined k-of-N architecture
%
function lambda_R = ...
    architecture_k_of_N( N, K, mu_l, lambda_s_x)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
TOL = 0.00000000000001;

Z_tab = lambda_r_to_p_z_x_marginal(lambda_s_x, mu_l); % compute 2*2 table. Z_tab(i,j) = Pr(X=i,X_R=j)
Z_tab = vec2mat(Z_tab,2);

mu =  sum(binopdf(K:N,N,mu_l)); % 1-binocdf(K-1,N,mu_l); % need at least K of N to be 'ON'

lambda_R = 0; % compute joint disease probability for two relatives
for j=K:N % number of 'on' liabilities in person
    for k=K:N % number of 'on' liabilities in relative
        for l=max(0,k+j-N):min(j,k) % number of shared 'on' liabilities
            lambda_R = lambda_R + Z_tab(2,2)^l * Z_tab(1,2)^(k+j-2*l) * Z_tab(1,1)^(N-k-j+l) * ...
                multinomial(N, [l, k-l, j-l, N-k-j+l]); % k+j-2*l
        end
    end
end
lambda_R = lambda_R ./ mu^2; % normalize by prevalence squared

