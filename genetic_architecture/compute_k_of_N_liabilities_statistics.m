% Compute all statistics for a special architecture which requires K of N liabilities.
% Computation is analytic.
%
% Input:
% N - number of loci
% K - minimal number of 'on' loci needed
% mu_l - 'prevalence' of each liability
% h_x - heritability (%genetic) of each liability
% h_shared_env - shared environment (%) in each liability. This is the fraction of TOTAL variance due to shared environemnt
% max_family_degree - how many generations away to compute 
% h_pop_str - which population estimator to use (default is 'ACE')
% 
% Output:
% lambda_R - risk for familial relatives of an individual
% stat_struct - structure with all statistical properties of architecture
%
function [lambda_R stat_struct] = ...
    compute_k_of_N_liabilities_statistics(N, K, mu_l, h_x, h_shared_env, ...
    max_family_degree, h_pop_str)

TOL = 0.00000000000001;
if(~exist('max_family_degree', 'var') || isempty(max_family_degree))
    max_family_degree = 6;
end
if(~exist('h_pop_str', 'var') || isempty(h_pop_str))
     h_pop_str = {'ACE'};
end
if(max_family_degree <0) % compute only the last lambda_R
    max_family_degree=-max_family_degree;
    keep_last_lambda_R=1;
else
    keep_last_lambda_R=0;
end

mu =  sum(binopdf(K:N,N,mu_l)); % 1-binocdf(K-1,N,mu_l); % need at least K of N to be 'ON'
x_mu_l = norminv(1-mu_l);
V = mu.*(1-mu);
lambda_R = zeros(max_family_degree,1); k_R_vec = [1 0.5 0.25 0.125 0.0625 0.03125];
lambda_type = {{'MZ-twins'}, {'DZ-twins', 'sibs', 'parent-offspring'}, ...
    {'grandparent-grandchild', 'half-sibs', 'uncle/ant'}, ...
    {'1st-cousins', 'graet-grandparent-great-grand-child'}, ...
    {'1st-cousing-once-removed', 'great-great-grandfother'}, {'2nd-cousings'}};

if(length(h_shared_env) == 1)% New! Allow different shared environment
    h_shared_env = repmat(h_shared_env, max(2,max_family_degree), 1);
end
for i=1:max(2,max_family_degree) %length(k_R_vec)
    lambda_R(i) = compute_lambda_R_LP_internal(N, K, h_x, h_shared_env(i), k_R_vec(i), mu_l, mu);    
end
lambda_MZ = lambda_R(1);
if(max_family_degree > 1)
    lambda_s = lambda_R(2);
else
    lambda_s = -1;  % temp indication that it's not computed
end

if(nargout > 1) % do heavy computations for all statistics
    corr_xz_second_term = -nchoosek(N-1, K-1) * mu_l^(K-1) * (1-mu_l)^(N-K); % no summing, just one term
    h = (N/(mu*(1-mu))) * (corr_xz_second_term * ...
        quadl(@(x)narrow_liability_conversion_fun(x, x_mu_l, sqrt(h_x)), ...
        -10*max(h_x, 1-h_x), 10*max(h_x, 1-h_x), TOL))^2; % compute narrow sense heritability on disease scale
    [h_liab_loci] = heritability_scale_change(h, 'liability', mu); % 'h_liab_loci_exact', h_liab_loci_exact,
    
    [h_liab_twins r_MZ r_DZ] = twin_concordance_to_heritability(lambda_MZ, lambda_s, mu, 'ACE');
    h_liab_pop = zeros(length(h_pop_str),1); 
    for j=1:length(h_pop_str)
       h_liab_pop(j) =  twin_concordance_to_heritability(lambda_MZ, lambda_s, mu, h_pop_str{j});
    end
    pi_liab_phantom = (h_liab_pop - h_liab_loci) ./ h_liab_pop;  % copmpute phantom heritability 
    
    h_liab_from_lambda_s = familial_risk_to_heritability(lambda_s, 'liability', mu, 0.5);
    h_liab_from_lambda_MZ = familial_risk_to_heritability(lambda_MZ, 'liability', mu, 1);
    h_liab_overestimation = (min(1,h_liab_twins) - h_liab_loci) / h_liab_loci; % this is in fraction of heritability overestimated out of true
    h_liab_unexplained_gap = (h_liab_twins - h_liab_loci) / h_liab_twins; % this is in fractions of heritability overestimated out of observed
    h_liab_unexplained_gap_from_MZ = (h_liab_from_lambda_MZ - h_liab_loci) / h_liab_from_lambda_MZ; % this is in fractions of heritability overestimated out of observed
    h_phantom_vec = h_liab_unexplained_gap; % New: call it in another name 
    

    h01_loci2 = N*h_x*normpdf(x_mu_l)^2 / (mu*(1-mu)^((2-N)/N)); % New! Compute heritability on the disease scale!! (works for K=1).
    
    
    H01 = mz_twin_risk_to_heritability(lambda_MZ, mu); % broad sense heritability on the disease scale
    h01_loci = h; % narrow sense heritability explained by loci on diseae scale
    r_MZ01 = mu .* (lambda_R(1) - 1) ./ (1-mu); % compute familial correlation on disease scale
    r_DZ01 = mu .* (lambda_R(2) - 1) ./ (1-mu);
    h01_pop = 2.*(r_MZ01-r_DZ01);
    h01_phantom_vec = 1 - h01_loci ./ h01_pop;
    
    stat_struct = struct('mu', mu, 'V', V, ...
        'lambda_s', lambda_s, 'lambda_MZ', lambda_MZ, ...
        'lambda_R', lambda_R, ... % 'lambda_type', lambda_type, ...
        'r_MZ', r_MZ, 'r_DZ', r_DZ, ...
        'h', h,  'h_liab_loci', h_liab_loci, ...
        'h_liab_twins', h_liab_twins, ... % 'h_liab_loci_exact', h_liab_loci_exact,
        'h_liab_from_lambda_s', h_liab_from_lambda_s, 'h_liab_from_lambda_MZ', h_liab_from_lambda_MZ, ...
        'h_liab_overestimation', h_liab_overestimation, 'h_liab_unexplained_gap', h_liab_unexplained_gap, ...
        'h_phantom_vec', h_phantom_vec, ...
        'h_liab_unexplained_gap_from_MZ', h_liab_unexplained_gap_from_MZ, ...
        'H01', H01, 'h01_loci', h01_loci, 'r_MZ01', r_MZ01, 'r_DZ01', r_DZ01, ...
        'h01_loci2', h01_loci2, ... % New! for debug
        'h01_pop', h01_pop, 'h01_phantom_vec', h01_phantom_vec, ...
        'h_liab_pop', h_liab_pop, 'pi_liab_phantom', pi_liab_phantom);
    
    stat_struct.lambda_type = lambda_type;
end

if(keep_last_lambda_R)
    lambda_R = lambda_R(end);
end

if(max_family_degree < length(lambda_R))
    lambda_R = lambda_R(1:max_family_degree);
end


