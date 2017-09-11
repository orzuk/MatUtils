% Compute marginal and pairwise effect size of two loci in the MLT model
%
% Input:
% h_x - heritability explained in EACH liability by each locus
% freq - MAF of each locus
% K - number of liabilites required to exceed the threshold
% N - number of total liabilities
% mu - prevalence (of DISEASE!)
%
% Output:
% p_x_y - joint probability of genotype and ONE LIABILITY
% p_x_z  - joint probability of genotype and PHENOTYPE
% p_x_x_z - joint porbability of TWO genotypes in DIFFERENT liabilities and PHENOTYPE
% disease_grr - genetic relative risk (for DISEASE) of genotype
%
function [p_x_y p_x_z p_x_x_z disease_grr] = heritability_to_p_z_x_MLT(h_x, freq, K, N, mu)

options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance
mu_l = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu)), 0, 1, options); % find mu_l that keeps the prevalence
GRR_liability = heritability_to_genetic_relative_risk(...
    h_x, 'liability', freq, mu_l); % genetic relative risk of one locus. Divide heritability by two to account for two alleles

p_x_y = genetic_relative_risk_to_p_z_x_marginal(freq, GRR_liability, mu_l); % joint probability of locus and liability
p_total_tab = compute_k_of_N_joint_tab(N, K, mu_l); % Pr. one liability and three liabilities on)

p_x_y = vec2mat(p_x_y', 2); % Compute joint distribution table of 0/1 alleles and phenotype
p_x_z = zeros(2); % joint probability of locus and disease
p_x_z(2,2) = p_x_y(1,2)*p_total_tab(2,1) / sum(p_total_tab(:,1)) + ...
    p_x_y(2,2)*p_total_tab(2,2) / sum(p_total_tab(:,2));
p_x_z(2,1) = freq - p_x_z(2,2); % locus is 'on'
p_x_z(1,2) = mu - p_x_z(2,2);
p_x_z(1,1) = 1 - p_x_z(2,2) - p_x_z(1,2) - p_x_z(2,1);
[~, disease_grr] = ...
    p_z_x_marginal_to_genetic_relative_risk(vec2row(mat2vec(p_x_z')));

p_y_given_x_one = p_x_y(2,2) / sum(p_x_y(:,2)); % probability of liability being 'ON' given the locus
p_y_given_x_zero = p_x_y(2,1) / sum(p_x_y(:,1));
%p_z_given_x_one = p_x_z(2,2) / sum(p_x_z(2,:));
%p_z_given_x_zero = p_x_z(1,2) / sum(p_x_z(1,:));


switch N % N=1 recieves a special treatment (both loci in the same liability)
    case 1
%         GRR_liability_both_loci = heritability_to_genetic_relative_risk(...
%             h_x, 'liability', freq, mu_l); % genetic relative risk of one locus. Divide heritability by two to account for two alleles
%         p_x_y_both_loci = genetic_relative_risk_to_p_z_x_marginal(freq, GRR_liability_both_loci, mu_l); % joint probability of locus and liability
        x_mu = norminv(1-mu); % set threshold 
        beta = sqrt(h_x ./ (freq.*(1-freq)));
        p_z_given_x1_x2(2,2) = 1-normcdf( (x_mu-2*beta*(1-freq)) ./ sqrt(1-2*h_x) ); % both genotypes are one
        p_z_given_x1_x2(2,1) = 1-normcdf( (x_mu-beta*(1-freq)+beta*freq) ./ sqrt(1-2*h_x) ); % one is zero, one is one
        p_z_given_x1_x2(1,2) = p_z_given_x1_x2(2,1); % 
        p_z_given_x1_x2(1,1) = 1-normcdf( (x_mu+2*beta*freq) ./ sqrt(1-2*h_x) ); % both genotypes are zero
        
    otherwise
        switch K % Probability of disease given each locus. Assume N>1 !!!
            case 1 % 'OR' gate
                p_z_given_x1_x2(2,2) = 1 - (1-mu_l).^(N-2).*(1-p_y_given_x_one).^2; % both genotypes are one
                p_z_given_x1_x2(2,1) = 1 - (1-mu_l).^(N-2).*(1-p_y_given_x_one).*(1-p_y_given_x_zero); % one is zero, one is one
                p_z_given_x1_x2(1,2) = p_z_given_x1_x2(2,1);
                p_z_given_x1_x2(1,1) = 1 - (1-mu_l).^(N-2).*(1-p_y_given_x_zero).^2; % both are zero
                
            case N % 'AND' gate
                p_z_given_x1_x2(2,2) = mu_l.^(N-2).*p_y_given_x_one.^2; % both genotypes are one
                p_z_given_x1_x2(2,1) = mu_l.^(N-2).*p_y_given_x_one.*p_y_given_x_zero; % one is zero, one is one
                p_z_given_x1_x2(1,2) = p_z_given_x1_x2(2,1);
                p_z_given_x1_x2(1,1) = mu_l.^(N-2).*p_y_given_x_zero.^2; % both are zero
        end % switch K
end % switch N

p_x_x = [1-freq freq]' * [1-freq freq]; % joint probability of x_1 and x_2 (under independence - no linkage)
p_x_x_z_one = p_z_given_x1_x2 .* p_x_x; % Pr(x_1,x_2,z=1)
p_x_x_z_zero = (1-p_z_given_x1_x2) .* p_x_x; % Pr(x_1,x_2,z=0)
p_x_x_z = mat2vec([mat2vec(p_x_x_z_zero) mat2vec(p_x_x_z_one)]')';


%p_z_given_x1_x2_null(2,2) = p_z_given_x_one^2/mu; % null expectation for x_1 and x_2
%p_z_given_x1_x2_null(2,1) = p_z_given_x_one*p_z_given_x_zero/mu;
%p_z_given_x1_x2_null(1,2) = p_z_given_x1_x2_null(2,1);
%p_z_given_x1_x2_null(1,1) = p_z_given_x_zero^2/mu;
%p_x_x_z_one_null = p_z_given_x1_x2_null .* p_x_x; % prob x,x and z=1
%p_x_x_z_zero_null = (1-p_z_given_x1_x2_null) .* p_x_x; % prob x,x and z=1
% Format P(x_1, x_2, z):
% [ p_000 p_001 p_010 p_011 p_100 p_101 p_110 p_111]

