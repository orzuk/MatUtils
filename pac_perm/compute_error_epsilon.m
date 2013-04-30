% A lot of calculations, to compute the confidence level delta we have for
% a given nsamples and error epsilon
% a certain overlap 1-\eps with confidence 1-\delta
% Note: delta can be a vector
function epsilon_obtained = compute_error_epsilon(rand_flag, one_side_flag, true_corr_flag, sigma, alpha, delta, nsamples, Ngenes, miu, prior)
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
TOL = 0.000001;
I = sqrt(-1);

TWO =2;

% Compute f_star
[f_star f_sigma x_alpha ] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, sigma, nsamples, alpha,miu,prior);  % We use here the approximation of Fisher
f_sigma = f_sigma / sqrt(Ngenes);
epsilon_obtained =  1 - (norminv(delta) * f_sigma + f_star); % prob. to be below the threshold  ?????

