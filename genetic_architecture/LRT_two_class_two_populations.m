% Compute likelihood-ratio-test statistic between model 
% H0: both populations have same s value
% H1: populations have different s values 
% Statistic used: Likeihood-ratio-test statistic: 
% 
function LRT = LRT_two_class_two_populations(D1, k_vec1, n_vec1, LL1, D2, k_vec2, n_vec2, LL2)

% Assume that the parameters for a single population were already fitted 
LRT = 2*(LL1 + LL2); 


% Now fit together the data from both populations. Note: this is tough because for each one we've got 
% a different demography !!!! can we maximize together? depends how maximizing log-likelihood works. 
(Can be too slow to do grid-search) 
s_vec = logspace(???); 

LL_joints = ???; 


LRT = LRT - 2*LL_joint; 
