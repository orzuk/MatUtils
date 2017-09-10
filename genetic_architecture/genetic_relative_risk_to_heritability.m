% Compute heritability assuming an additive/multiplicative (?) model
%
% Input:
% f_vec - risk allele frequencies
% GRR_marginal - genetic relative risk for each SNP
% mu - disease prevalence
% model - (optional) multiplicative or additive (or logistic?)
%
% Output:
% lambda_s_vec - disease risk for a sibling of a diseased person
% lambda_s_add - disease risk for a sibling for additive model
% lambda_mz_add - disease risk for monozigotic twin for additive model
% h_add - fraction of heritability explained by an additive model
% V_add - additive variance of trait
% lambda_s_mult - disease risk for a sibling for multiplicative model
% lambda_mz_mult - disease risk for monozigotic twin for multiplicative model
% h_mult - fraction of heritability explained by a multiplicative model
% V_mult - multiplicative variance of trait
% sib_freq_mat - matrix of joint disease freq. for person and sibling
% h_liab - heritability on liability scale
% h_liab_vec - heritability on liability scale of each locus 
% 
function [lambda_s_vec ...
    lambda_s_add lambda_mz_add h_add V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab h_liab_vec] = ...
    genetic_relative_risk_to_heritability(f_vec, GRR_marginal, mu, model, varargin)

AssignGeneralConstants;
f_vec = vec2column(f_vec);
sib_freq_mat = [];
n = length(f_vec);
lambda_s_vec = zeros(n,1); % For lambda_s we need only the odds ratios

if(mod(n,2) == 0) % This is just for checking diploid computations. Results are not really used. 
    for i=1:2:n-1 % loop on SNPs. Each SNP represented by two genotypes
        i_is = i; % compute lambda_s - method from Eliana
        q = f_vec(i); % get MAF
        T = [1-q+q^2/4, q-q^2/2, q^2/4; ...
            1/2-3*q/4+q^2/4, 1/2+q/2-q^2/2, q/4+q^2/4; ...
            (1-q)^2/4, 1/2-q^2/2, 1/4+q/2+q^2/4]; % transfer matrix
        D = diag( [(1-q)^2, 2*q*(1-q), q^2] ); % genotype frequencies matrix
        v = [1 GRR_marginal(i) GRR_marginal(i)*GRR_marginal(i+1)]; % assume multiplicative model
        lambda_s_vec(i) = v * D * T * v'; %  / mu^2; % don't divide by mu^2 - this is accounted for in v
        sib_freq_mat = D * T;
    end
    lambda_s_vec = lambda_s_vec ./ ( (1 - f_vec).^2 + ...
        2 .* f_vec .* (1-f_vec) .* GRR_marginal + ...
        f_vec .^ 2 .* GRR_marginal .^ 2).^2; % correct for prevalence (assume all are the same)
    lambda_s2 = prod(lambda_s_vec(1:2:end-1));
end

% Try primitive version: work on each allele seperately
for i=1:n % loop on SNPs. Each SNP represented by two genotypes
    i_is = i; % compute lambda_s - method from Eliana
    q = f_vec(i); % get MAF
    T = [1 - q/2, q/2; (1-q)/2, (1+q)/2]; % conditional genotype probability matrix T(i,j) = Pr(x_s = j | x=i)
    D = diag( [1-q, q] ); % allele frequencies matrix
    v = [1 GRR_marginal(i)]; % assume multiplicative model
    lambda_s_vec(i) = v * D * T * v'; %  / mu^2; % don't divide by mu^2 - this is accounted for in v
end
lambda_s_vec = lambda_s_vec ./ (1 - f_vec + f_vec .* GRR_marginal).^2; % correct for prevalence
lambda_s = prod(lambda_s_vec); %  ./ prod(1 - f_vec + f_vec .* odds_ratio_marginal).^2;
lambda_s_mult = lambda_s; % same value for different named variable
if(mod(n,2) == 0) % check diploid computations
    both_lambdas_should_be_equal = (lambda_s2 - lambda_s) /  max(lambda_s2, lambda_s); % compute relative difference (lambda could be huge)
    if(abs(both_lambdas_should_be_equal) > 0.000000001) % allow a little greater tollerance epsilon)
        both_lambdas_should_be_equal = (lambda_s2 - lambda_s) /  max(lambda_s2, lambda_s) % compute relative difference (lambda could be huge)
        error('why different lambdas???');
    end
end

family_tree = generate_family_tree(3); K = kinship_coefficient(family_tree); K_s = K(end-1,end); % get kinship coefficient
T_11 = f_vec .* (K_s + (1-K_s) .* f_vec); % get transfer matrix prob.: Pr(x_i=1, x_S,i = 1)
alpha_vec = mu .* (GRR_marginal - 1) ./ (1 - f_vec + f_vec .* GRR_marginal); % alpha_vec = GRR_marginal;
beta = mu - sum(alpha_vec .* f_vec);
lambda_s_add = (alpha_vec .* f_vec) * (alpha_vec .* f_vec)'; % take interaction term when assuming additive model
lambda_s_add = sum(sum(lambda_s_add - diag(diag(lambda_s_add)))) + ...
    beta^2 + 2*beta * sum(alpha_vec .* f_vec) + sum(alpha_vec .^ 2 .* T_11); % compute via additive model
lambda_s_add = lambda_s_add / mu^2; % normalize by indepdnent relative risk

lambda_mz_mult = prod(1 - f_vec + f_vec .* GRR_marginal.^2) ./ ...
    prod(1 - f_vec + f_vec .* GRR_marginal).^2; % Compute risk for monzygotic twins (assuming multiplicative model)

V_mult = mu*(1-mu); % for additive heritability we need also prevalence
V_add = f_vec .* (1-f_vec) .* alpha_vec.^2; % new: marginal additive variance in additive model
h_add = sum(V_add) / V_mult;

V_env = mu  - mu^2 * lambda_mz_mult; h_mult = 1 - V_env / V_mult;
lambda_mz_add = heritability_to_mz_twin_risk(h_add, mu); % this includes only additive variance! 

h_liab = heritability_scale_change(h_add, 'liability', mu); %h_add * V_mult / z_score^2; % Compute also liability threshold model heritability
h_liab_vec = heritability_scale_change(V_add ./ V_mult, 'liability', mu); % liability-scale heritability of each variant 
