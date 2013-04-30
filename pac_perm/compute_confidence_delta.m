% A lot of calculations, to compute the confidence level delta we have for
% a given nsamples and error epsilon
% a certain overlap 1-\eps with confidence 1-\delta
% Note: epsilon can be a vector
function delta_obtained = compute_confidence_delta(rand_flag, one_side_flag, true_corr_flag, sigma, alpha, epsilon, nsamples, Ngenes, miu, prior)
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
TOL = 0.000001;
I = sqrt(-1);
% miu = []; prior = [];

TWO =2;

% Compute f_star
[f_star f_sigma x_alpha ] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, sigma, nsamples, alpha,miu,prior)  % We use here the approximation of Fisher
f_sigma = f_sigma / sqrt(Ngenes);
delta_obtained =  normcdf((1-epsilon-f_star)/f_sigma); % prob. to be below the threshold 

