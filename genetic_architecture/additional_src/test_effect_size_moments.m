% Compute the distribution of effect size estimated from a set of samples
function [effect_mu effect_std effect_mu_min effect_mu_max] = ...
    test_effect_size_moments(n_cases, n_controls, effect_size, mu, f, effect_type, alpha, iters)

if(~exist('iters', 'var') || isempty(iters)) % set confidence interval
    iters = 1000;
end

p_z_x = genetic_relative_risk_to_p_z_x_marginal(f, effect_size, mu);
p_z_x = pop_prob_to_case_control_prob(p_z_x, n_cases, n_controls);

n_samples = n_cases + n_controls;
contigency_table_vec = mnrnd(n_samples, p_z_x, iters); % randomize table
contigency_table_vec = normalize(contigency_table_vec, 2)   

effect_dist = p_z_x_marginal_to_genetic_relative_risk(contigency_table_vec); 

figure; hist(effect_dist, 500); % plot histogram of estimated effect sizes
effect_mu = mean(effect_dist);
effect_std = std(effect_dist);

if(~exist('alpha', 'var') || isempty(alpha)) % set confidence interval
    alpha = 0.05;
end
num_stds = norminv(1-alpha/2); % set # of st.d.s. to take
effect_mu_min = effect_mu - num_stds * effect_std;
effect_mu_max = effect_mu + num_stds * effect_std;

