% Compute joint probability distribution of one variable x_i and the entire
% output of a binomial distribution: z = 1_{\sum_i x_i >= k}
%
% Input:
% N - number of binary r.v.s.
% k - how many binary r.v.s. need to be 'turned on' (at least) to have z=1
% p - probability of each variable being on
%
% Output:
% jount_tab - a 2x2 table with probabilties for a genotype and phenotype Pr(x,z)
%
function joint_tab = compute_k_of_N_joint_tab(N, K, p)

joint_tab = zeros(2,2);

if(K == 1) % more accurate computation for low prevalence
    mu = sum(binopdf(K:N,N,p));
else
    mu = 1-binocdf(K-1,N,p); % probability at least K are 'on'
end
joint_tab(2,2) = p* (1-binocdf(K-2,N-1,p)); % both variables are set to one
joint_tab(1,2) = p-joint_tab(2,2); % disease=0, genotype=1
joint_tab(2,1) = mu-joint_tab(2,2); % disease=1, genotype=0
joint_tab(1,1) = 1-joint_tab(1,2)-joint_tab(2,1)-joint_tab(2,2); 
