% A  script for computing overlap for Heike's data 
AssignAllGlobalConstants;

% Choose the q probability distribution of the TRUE corrleations
rand_flag = GAUSSIAN;           % take q to be Gaussian
one_side_flag = ONE_SIDE;       % Choose if to take both top and bottom genes or just top
true_corr_flag  = TWO_SAMPLED;  % Choose if to compute overlap between true and sampled or between two sampled
const_alpha_flag = 0;           % Flag saying if we take alpha to be const or take alpha*N ~ n

% Here do only maple integral to show a pic.
res = 0.005; max_sigma = 4; alpha_res = 0.01;
sigma_vec = [res:res:max_sigma];

% Here there's only one alpha!
N_g = 34842; 
top_genes=100;
n_samples=65;
alpha = top_genes/N_g;          %alpha_vec = [alpha_res:alpha_res:0.5];
dist_std = 0.1198959;           % The std. of the original distribution q. If it's not 1 then everything re-scales.
                                % its sqrt(sigma^2)
miu = []; prior = [];

%%%%% Maple Integral 
%%%% int( (1/sqrt(2*Pi)) *  exp(-c*c/2) * int( (1/(sigma*sqrt(2*Pi))) * exp(-x*x/(2*sigma*sigma)),x=x_alpha-c..infinity), c=c_alpha..infinity);
C_alpha = norminv(1-0.5*alpha); % Get the C_alpha vector


% Here define the vector with no of samples
n_res = 5; max_n = 1000;
samples_vec=[n_res:n_res:max_n];
f_mean_vec = zeros(1,length(samples_vec));
f_std_vec = zeros(1,length(samples_vec));
% f_mean_vec = zeros(1,length(sigma_vec));
% f_std_vec = zeros(1,length(sigma_vec));
% f_vec = [res:res:1-res];
% nsamples_vec = 1./sigma_vec.^2 + 3;

for i = 1:length(samples_vec)
  [f_mean_vec(i), f_std_vec(i), x_alpha] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, dist_std, samples_vec(i), alpha, miu, prior);
end
figure; hold on; plot(samples_vec, f_mean_vec); 
