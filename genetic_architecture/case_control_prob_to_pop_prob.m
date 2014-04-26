% Transfer a joint case-control distribution to population distribution 
% (Pr(Z=1)= prevalence, frac. of cases)
%
% Input:
% p_vec - joint probability of one or two genotypes and phenotype
% prevalence - fraction of cases in population 
%
% Output:
% q_vec - adjusted porbability vector in the population 
%
function q_vec = case_control_prob_to_pop_prob(p_vec, prevalence)

q_vec = pop_prob_to_case_control_prob(p_vec, prevalence, 1-prevalence); % set proportions of cases and controls

