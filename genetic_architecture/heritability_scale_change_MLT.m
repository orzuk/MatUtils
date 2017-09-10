% Change heritability computed based on one LT model to heritability at each liability
% assuming MLT model
% 
% Input: 
% h_x - heritability explained by a locus for the trait, assuming LT model
% K - number of liabilities required to exceed the threshold
% N  - total number of liabilities
% mu - prevalence
% output_scale - conversion direction MLT->LT (default) or LT->MLT
% 
% Output: 
% h_liab_out - converted heritability on LT (default) or MLT scale
% 
function h_liab_out = heritability_scale_change_MLT(h_x, K, N, mu, output_scale)

TOL = 0.00000000000000001;
options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance

%h = heritability_scale_change(h_x, 'binary', mu); % compute narrow-sense heritability
mu_l = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu)), 0, 1, options); % find mu_l that keeps the prevalence
x_mu_l = norminv(1-mu_l); % compute threshold 

if(~exist('output_scale', 'var') || isempty(output_scale))
    output_scale = 'LT'; 
end
switch output_scale
    case 'LT' % compute the single LT heritability given the MLT heritability 
        corr_xz_second_term = -nchoosek(N-1, K-1) * mu_l^(K-1) * (1-mu_l)^(N-K); % no summing, just one term
        h_binary = (N/(mu*(1-mu))) * (corr_xz_second_term * ...
            quadl(@(x)narrow_liability_conversion_fun(x, x_mu_l, sqrt(h_x)), ...
            -10*max(h_x, 1-h_x), 10*max(h_x, 1-h_x), TOL))^2; % compute narrow sense heritability on disease scale
        h_liab_out = heritability_scale_change(h_binary, 'liability', mu);        
    case 'MLT' % Perform reveresed binary search. Takes a long time 
        options = optimset('tolx', TOL); % increase optimization tolerance         
        h_liab_out = fminbnd(@(x) abs(heritability_scale_change_MLT(x, K, N, mu, 'LT')-h_x), ...
            0,1, options);        
end
