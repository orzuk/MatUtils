% Determine for each MAF what is the minimal effect size discovered in a GWAS
%
% Input:
% x_vec - MAF
% mu - prevalence (used only for binary traits)
% num_cases - # of cases/total population in study
% num_controls - # of controls in study
% alpha - pvalue cutoff
% test_type - type of statistical test
% test_stat - test statistic used
%
% Output:
% critical_beta_vec - critical beta above which effect is declared as significant
% critical_beta_inv_vec - critical beta below which effect is declared as significant
%
function  [critical_beta_vec critical_beta_inv_vec] = compute_discovery_boundary(x_vec, mu, ...
    num_cases, num_controls, alpha, test_type, test_stat) % Determine discovery boundry

x_vec = vec2column(x_vec); % TEMP!!
AssignGeneralConstants;
num_loci = length(x_vec);

max_beta_vec = vec2column((1-epsilon) ./ sqrt(2 .* x_vec .* (1-x_vec)));
max_beta_vec = max(max_beta_vec, 1); % allow at least beta=1
if(~isempty(strfind(test_stat, 'QTL'))) % QTL
    trait_type = 'QTL';
    min_beta_vec = zeros(num_loci,1);
else % binary. This is actually GRR
    trait_type = 'binary';
    min_beta_vec = ones(num_loci,1);
    max_beta_vec = vec2column(0.5./ sqrt(2 .* x_vec .* (1-x_vec))); % take much lower values to avoid numerical problems with high heritability
    %    max_beta_vec = max(100, beta_to_genetic_relative_risk(max_beta_vec, x_vec, mu)); % guarantee at least 100 GRR
    max_beta_vec(:) = 1000; % quick and dirty
end
while(max(max_beta_vec - min_beta_vec) > sqrt(epsilon)) % perform binary search
    mid_beta_vec = (max_beta_vec + min_beta_vec) ./ 2;
    if(~isempty(strfind(test_stat, 'QTL')))
        p_vec = QTL_params_to_p_mat(x_vec, mid_beta_vec, 1);
    else % binary traits
        p_vec = genetic_relative_risk_to_p_z_x_marginal(x_vec, mid_beta_vec, mu);
    end
    cur_power_vec = compute_association_power(p_vec, ...
        num_cases, num_controls, alpha, 1, test_type, test_stat,[],1);
    high_inds = find(cur_power_vec == 1); low_inds = find(cur_power_vec == 0);
    min_beta_vec(low_inds) = mid_beta_vec(low_inds);
    max_beta_vec(high_inds) = mid_beta_vec(high_inds);
    
    %    .* (1-cur_power_vec) + ...
    %        max_beta_vec .* cur_power_vec);
end

critical_beta_vec = mid_beta_vec;
%min_cur_power_vec = compute_association_power(x_vec, ...
%    num_cases, num_controls, alpha, test_type, test_stat)

if(nargout > 1)
    switch trait_type
        case 'QTL'
            critical_beta_inv_vec = -critical_beta_vec; % simply take negative
        case 'binary'
            critical_beta_inv_vec = 1 ./ compute_discovery_boundary(1-x_vec, mu, ...
                num_cases, num_controls, alpha, test_type, test_stat); % 1/GRR of 1-MAF
    end
end




