% Script for testing computation of constraints on genetic architecutres from first few moments

N = 10000;
% s_vec =  -logspace(-6,-1,100);
s_vec = linspace(-0.1, -0.000001, 500); 
beta_vec = -1:0.01:1; %-s_vec;
coupling_str = 'uniform'; %  'linear';


coupling_str_vec = {'uniform', 'linear', 'sigmoid', 'eyre-walker'};
coupling_params_vec = {0, 0.1,  30}; % {'uniform', 'linear'};
coupling_params_vec{4}.tau = 1; coupling_params_vec{4}.sigma = 1; coupling_params_vec{4}.theta = 1; coupling_params_vec{4}.k = 1;
num_couplings = length(coupling_str_vec); 
mean_phenotype_change = zeros(num_couplings, 1); 
for i=1:num_couplings
    mean_phenotype_change(i) = compute_architecture_constraints(coupling_str_vec{i}, coupling_params_vec{i}, N, s_vec, beta_vec);
end
    
%mean_square_phenotype_change ...
%    skewness kurtosis heritability] = compute_architecture_constraints(coupling_str, N, s_vec, beta_vec)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_theta_vec = logspace(-4, 4, 1000); h=1/2;
beta_vec = sqrt(s_theta_vec * h/2); 
Kurt = 0.2; % example 
max_beta = sqrt(Kurt / h);
[~, max_ind] = find(beta_vec >= max_beta, 1); 
figure; 
loglog(s_theta_vec, beta_vec, ':', 'linewidth', 2); hold on; 
loglog(s_theta_vec, repmat(max_beta, length(s_theta_vec), 1), 'k--', 'linewidth', 2); 
loglog(s_theta_vec(1:max_ind), beta_vec(1:max_ind), 'linewidth', 3); hold on; 
%loglog(s_theta_vec, -beta_vec, 'linewidth', 3); 
xlabel('-S/\theta', 'fontsize', 16); ylabel('|\beta|', 'fontsize', 16); 
my_saveas(gcf, 'beta_vs_S', {'pdf', 'jpg', 'eps'}); 




figure; 
plot(s_theta_vec, repmat(max_beta, length(s_theta_vec), 1), 'k--'); 
plot(s_theta_vec, repmat(-max_beta, length(s_theta_vec), 1), 'k--'); 

plot(s_theta_vec, beta_vec, 'linewidth', 3); hold on; 
plot(s_theta_vec, -beta_vec, 'linewidth', 3); 

xlabel('-S/\theta', 'fontsize', 16); ylabel('\beta', 'fontsize', 16); 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%