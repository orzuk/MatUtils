% Compute chi-square non-centrality parameter for a locus (indicative of power)
% Problem: this depends on the effect size and RAF. There's no 'universal' power measure
% Note: The 'ott' option matches Purcell's genetic power calculator
% version: 'Case-control statistics: allelic 1 df test (B versus b)'
%
% Input:
% f_vec - risk allele frequency
% grr_vec - effect size (genetic relative risk)
% num_cases - # of cases in study
% num_controls - # of controls in study
% mu - prevalence (optional) should be used but can be neglected when small
% ncp_type - how to compute ncp parameter (several approximations are available)
% 
% Output:
% non_centrality_parameter - chi-square non-centrality. Marks deviation from null-hypothesis and indicates power
%
function non_centrality_parameter = genetic_relative_risk_to_non_centrality_parameter( ...
    f_vec, grr_vec, num_cases, num_controls, mu, ncp_type)

num_loci = length(f_vec);
num_samples = num_cases+num_controls; 
theta = num_cases / num_samples;
if(~exist('ncp_type', 'var') || isempty(ncp_type))
    ncp_type = 'standard_chi_square'; % 'apples';
end
if(~exist('mu', 'var') || isempty(mu)) % can't use mu
    ncp_type = 'apples';
end
switch lower(ncp_type)
    case 'standard_chi_square'
        num_loci = length(f_vec);
        non_centrality_parameter = zeros(num_loci,1);
        p_z_effect_vec = genetic_relative_risk_to_p_z_x_marginal(f_vec, grr_vec, mu);
        p_z_effect_vec = pop_prob_to_case_control_prob(p_z_effect_vec, num_cases, num_controls); % adjust to case-control!
        for i=1:num_loci
            p_z_marginal_vec = vec2row(mat2vec(table_to_marginal_probs(vec2mat(p_z_effect_vec(i,:), 2))));
            non_centrality_parameter(i) = num_samples * ...
                sum((p_z_effect_vec(i,:) - p_z_marginal_vec).^2 ./ p_z_marginal_vec);  % here a direct computation of general chi-square table (assume cases=controls)
        end
    case 'apples'
        non_centrality_parameter = 2*num_samples .* f_vec .* (1-f_vec) .* theta .* (1-theta) .* (grr_vec-1).^2; % log(grr_vec.^2);
        if(exist('mu', 'var'))
            non_centrality_parameter = non_centrality_parameter ./ ...
                ((1-mu).^2 * (1+f_vec.*(grr_vec-1)));
        end
    case 'eliana'
        non_centrality_parameter = 2*num_samples .* f_vec .* (1-f_vec) .* theta .* (1-theta) .* log(grr_vec).^2; % log(grr_vec.^2);
        
        % Should use 'ott' below - it matches Purcell's genetic power calculatur
    case 'ott' % from http://linkage.rockefeller.edu/pawe3d/help/Linear-trend-test-ncp.html
        if(num_cases == 0) % no way we can see signal here!!! 
            non_centrality_parameter = 0;
            return;
        end
        p_z_effect_vec = genetic_relative_risk_to_p_z_x_marginal(f_vec, grr_vec, mu);
        p_z_effect_vec = pop_prob_to_case_control_prob(p_z_effect_vec, num_cases, num_controls); % adjust to case-control!
        %        p_z_effect_vec = mat2vec(expand_allele_table(vec2mat(p_z_effect_vec, 2)', 'trend'))';  % Need to check and make sure expansion is still valid here !!!!
        p_z_effect_vec = expand_allele_table(p_z_effect_vec, 'trend', '4XN'); % cases inds are 2 4 6
        p_given_case_vec = p_z_effect_vec(:,[2 4 6]) ./ repmat(sum(p_z_effect_vec(:,[2 4 6]),2), 1, 3); % prob. of genotype given that we picked a case
        p_given_control_vec = p_z_effect_vec(:,[1 3 5]) ./ repmat(sum(p_z_effect_vec(:,[1 3 5]),2), 1, 3); % prob. of genotype given that we picked a control
        non_centrality_parameter = num_cases*num_controls * ...
            sum(repmat([0 1 2], num_loci, 1) .* (p_given_case_vec - p_given_control_vec), 2).^2;
        non_centrality_parameter = non_centrality_parameter ./ ...
            ( sum(repmat([0 1 2].^2, num_loci, 1) .* ...
            (num_cases .* p_given_case_vec + num_controls .* p_given_control_vec), 2) - ...
            sum(repmat([0 1 2], num_loci, 1) .* ...
            (num_cases .* p_given_case_vec + num_controls .* p_given_control_vec), 2).^2./num_samples ); % normalize
end




% compute.power <- function(grr, freq, p.val, Ncase, Ncontrol){
% #First compute noncentrality parameter for the trend test
% N <- Ncase+Ncontrol;
% theta <- Ncase/N;
% nc.param <- 2*N*freq*(1-freq)*theta*(1-theta)*(log(grr)^2);
% k = 1 #1 degree of freedom
% z <- qchisq(1-p.val,1)
% power <- 1-pchisq(z,1,ncp=nc.param);
% return(power);
% }
%
% compute.power.QTL <- function(q.sq,N){
%     z <- qchisq(1-5e-8,1)
%     NCP <- N*q.sq/(1-q.sq);
%     power <- 1-pchisq(z,1,ncp=NCP);
%     return(power);
%     }
%
%     [Hide Quoted Text]