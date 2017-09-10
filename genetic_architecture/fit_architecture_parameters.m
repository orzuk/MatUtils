% Fit the best additive/multiplicative model for a certain architecture
% Models assume no interactions but some transofmation from 'risk' to
% disease - so they are essentially generalized additive models.
% The fitted model may give architecture with different mean (prevalence)
% than the true one - thus they are not neccessarily unbiased estimators.
%
% Input:
% params_struct - structure with architecture parameters
% architecture_str - type of architecture
% f_vec - frequencies of each allele
% fit_mode - what function to fit
% GRR_marginal - genetic relative risk for each allele
% compute_method_flag - how to generate inputs&outputs for architecture (currently only enumerate or sampling)
%
% Output:
% alpha - linear coefficients vec
% beta - affine constant
% structure representing the fitted architecture
%
function [alpha beta fitted_params_struct] = fit_architecture_parameters(params_struct, ...
    architecture_str, f_vec, fit_mode, GRR_marginal, compute_method_flag)

AssignStatsConstants(); 
if(~exist('compute_method_flag', 'var') || isempty(compute_method_flag))
    compute_method_flag = 'enumerate';
end
N = length(f_vec);

switch compute_method_flag
    case 'enumerate'
        if(N > 20)
           error('Cant enumerate! number of SNPs too large!'); 
        end        
        iters = 1; 
    case 'sampling'
        iters = 10000; % number of individuals to sample 
    case 'analytic'
        compute_method_flag = 'sampling'; % we don't know yet to fit analytically
        iters = 10000; 
end
[x_vec p_x_vec] = initilize_x_vec_constants(N, 0, f_vec, compute_method_flag, iters); % generate input and output vec
z = genetic_architecture(x_vec, architecture_str, params_struct, 1, BINARY); % sample zero-one genotypes
mu = sum(z .* p_x_vec);
% weights_vec = repmat(sqrt(p_x_vec), 1, N+1)
x_vec_reg = [ ones(size(x_vec,1),1) x_vec] .* repmat(sqrt(p_x_vec), 1, N+1); % normalize by sqrt for regression
switch compute_method_flag
    case 'sampling' % here all sampled vectors are weihted the same
        p_x_vec(:) = 1 ./ iters;
end

switch fit_mode
    case 'additive_naive'
        mu_marginal = mu .* vec2row(GRR_marginal) ./ (1 - f_vec + f_vec .* vec2row(GRR_marginal));
        alpha = mu_marginal - (mu - f_vec .* mu_marginal) ./ (1 - f_vec);
        alpha = [mu - sum(alpha .* f_vec) alpha];
        
    case 'additive'
        take_marginals = 0;
        if(take_marginals)
            mu_marginal = mu .* GRR_marginal ./ (1 - f_vec + f_vec .* GRR_marginal);
            alpha = mu_marginal - (mu - f_vec .* mu_marginal) ./ (1 - f_vec);
            alpha = [mu - sum(alpha .* f_vec) alpha];
        else
            z_vec_reg = z .* sqrt(p_x_vec);
            [alpha alpha_int alpha_residuals] = regress(z_vec_reg, x_vec_reg);
        end
        
    case 'multiplicative_naive'
        alpha = vec2row(log(GRR_marginal)); %  beta = log(min(z)); % take the minimal possible value
        alpha = [log(mu) - sum( log(1-f_vec + f_vec .* exp(alpha))) alpha];
    case 'multiplicative' % Here we simply take the marginal relative risk of each one (?)
        take_marginals = 0;
        if(take_marginals) % just take genetic relative risk of each allele separetely
            alpha = vec2row(log(GRR_marginal)); %  beta = log(min(z)); % take the minimal possible value
            alpha = [log(mu) - sum( log(1-f_vec + f_vec .* exp(alpha))) alpha];
        else % really fit a multiplicative model
            z_vec_reg = z .^ sqrt(p_x_vec); % multiply in log space
            %            alpha = glmfit(x_vec_reg(:,1:end-1), z_vec_reg, 'binomial', 'link', 'log');
            [alpha dev] = glmfit(x_vec, z, 'normal', 'link', 'log', 'weights', p_x_vec); % fit with weights for each observation
            %            [alpha dev] = glmfit(x_vec, z, 'normal', 'link', 'log'); % fit with weights for each observation
        end
    case 'logistic'
        alpha = glmfit(x_vec, z, 'binomial', 'link', 'logit', 'weights', p_x_vec); % fit with weights for each observation
        
    case {'liability', 'liability-threshold', 'probit'} % fit the liability threshold model often used in genetics
        alpha = glmfit(x_vec, z, 'binomial', 'link', 'probit', 'weights', p_x_vec); % fit with weights for each observation
    otherwise
        alpha = -999999999; fitted_params_struct = params_struct;
end
beta = alpha(1); alpha = alpha(2:end); % seperate constant term from others 

if(~exist('fitted_params_struct', 'var')) % generate a structure with fitted architecture 
    fitted_params_struct = my_rmfield(params_struct, ...
        {'gate_type', 'gate_params', 'circuit', 'circuit_gate_params'});
    fitted_params_struct.num_clauses = 1;
    fitted_params_struct.k_in_clause = N;
    fitted_params_struct.linear_coef_vec = alpha;
    fitted_params_struct.b = beta;
    fitted_params_struct.architecture_str = fit_mode;
    if(~isfield(fitted_params_struct, 'z_std'))
        fitted_params_struct.z_std = 0;
    end
end
