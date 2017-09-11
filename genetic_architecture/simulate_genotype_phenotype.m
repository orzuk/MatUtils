% Simulate genotypes for ONE SNP and phenotypes for a set of individuals
%
% Input:
% p_vec - matrix with genotype/phenotype frequencies
% num_cases - sample size (cases, or population for QTL)
% num_controls - sample size (controls)
% iters - number of output vectors
% trait_type - binary or QTL
% trait_mode - can be genotypes (2x2) or alleles (3x2) table (should enable also epistasis 4x2)
% full_flag - 0 (default) returns only summary statistics, 1 returns a vector of genotype and phenotype
%
% Output:
% contigency_table - table of genotype and phenotype frequencies. If full flag is on, output is 2 vectors: genotype and phenotype
% genotype_phenotype_prod - sum_i x_i y_i (needed as another sufficient statistic)
% genotype_sqr - sum_i x_i^2 (needed as another sufficient statistic)
%
function [contigency_table genotype_phenotype_prod genotype_sqr] = ...
    simulate_genotype_phenotype(p_vec, num_cases, num_controls, ...
    iters, trait_type, trait_mode, full_flag, varargin)

snp_type = 'diploid'; % default: diploid snps
if(~exist('iters', 'var') || isempty(iters))
    iters = 1;
end
if(~exist('trait_type', 'var') || isempty(trait_type)) % set default trait type
    trait_type = 'QTL';
end
if(~exist('trait_mode', 'var') || isempty(trait_mode))
    trait_mode = 'allele'; % default is two alleles
end
if(~exist('full_flag', 'var') || isempty(full_flag))% set default: just summary statistics
    full_flag = 0;
end

switch trait_type     % Binary trait or QTL
    case {'binary', 'Binary'}
        if(full_flag) % here actually simulated a vector of genotype/phenotype
            num_samples = num_cases + num_controls; % Temp! simulate cases and controls together !!! 
            [f_vec beta_vec mu] = p_mat_to_QTL_params(p_vec); % extract from encoding effect size and MAF
            N = p_vec(3); K = 1; % Temp. encoding 
            h_x = 2*f_vec*beta_vec.^2; % heritability of a locus (on liability scale)            
            mu_l = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu)), 0, 1); % find mu_l that keeps the prevalence
            x_mu = norminv(1-mu_l); 
            switch trait_mode
                case {'genotype', 'genotypes'}
                    % MISSING HERE!!
                    contigency_table = zeros(iters, num_samples,2); % For each iteration and each individual generate x1, z, z_expected
                case {'allele', 'alleles'} % here there's no meaning: it's liability
                    % MISSING HERE!!
                    contigency_table = zeros(iters, num_samples,3); % For each iteration and each individual generate x1, z, z_expected
                case {'epistasis', 'interaction', 'pairwise'} % Here Gaussian r.v.s.
                    contigency_table = zeros(iters, num_samples,4); % For each iteration and each individual generate x1, x2, z, z_expected
                    contigency_table(:,:,1) = randn(iters, num_samples); % x1 normalized .* sqrt(h_x); % x1
                    contigency_table(:,:,2) = randn(iters, num_samples); % x2 normalized .* sqrt(h_x); % x2
                    
                    if(N >= 2) % loci are in different liabilities
%                        contigency_table(:,:,3) = rand(iters, num_samples) > ... % z
                         contigency_table(:,:,4) =  1 - ((1-mu_l).^(N-2) .* ...
                            normcdf( (x_mu - contigency_table(:,:,1) .* sqrt(h_x)) ./ sqrt(1-h_x)) .* ...
                            normcdf( (x_mu - contigency_table(:,:,2) .* sqrt(h_x)) ./ sqrt(1-h_x)));  % z_expected
                    else % loci are in the same liability 
%                        contigency_table(:,:,3) = rand(iters, num_samples) > ... % z
                        contigency_table(:,:,4) = 1 - (normcdf( (x_mu - contigency_table(:,:,1) .* sqrt(h_x) - ...
                            contigency_table(:,:,2) .* sqrt(h_x)) ./ ...
                        sqrt(1-2*h_x))); % z_expected
                    end
                    contigency_table(:,:,3) = rand(iters, num_samples) < contigency_table(:,:,4);
            end
            genotype_phenotype_prod = []; genotype_sqr = []; % not needed since we've got the full data
%            var_z = contigency_table(:,:,3); % Compute Variance explained
                        
        else % here just simulate summary tables            
            switch trait_mode
                case {'genotype', 'genotypes'}
                    contigency_table = zeros(iters, 6);
                    contigency_table(:,[1 3 5]) = mnrnd(num_controls, p_vec([1 3 5])./sum(p_vec([1 3 5])), iters);
                    contigency_table(:,[2 4 6]) = mnrnd(num_cases, p_vec([2 4 6])./sum(p_vec([2 4 6])), iters); % randomize counts in cases and controls
                case {'allele', 'alleles'}
                    contigency_table = zeros(iters, 4);
                    contigency_table(:,[1 3]) = mnrnd(num_controls, p_vec([1 3])./sum(p_vec([1 3])), iters);
                    contigency_table(:,[2 4]) = mnrnd(num_cases, p_vec([2 4])./sum(p_vec([2 4])), iters); % randomize counts in cases and controls
                case {'epistasis', 'interaction', 'pairwise'}
                    contigency_table = zeros(iters, 8);
                    contigency_table(:,[1 3 5 7]) = mnrnd(num_controls, p_vec([1 3 5 7])./sum(p_vec([1 3 5 7])), iters);
                    contigency_table(:,[2 4 6 8]) = mnrnd(num_cases, p_vec([2 4 6 8])./sum(p_vec([2 4 6 8])), iters);
            end % switch trait mode
        end % if full
    case {'quantitative', 'QTL'} % phenotype is QTL
        [f_vec beta_vec] = p_mat_to_QTL_params(p_vec); % extract from encoding
        mog_mu = beta_vec; % 1; % beta_vec ./ sqrt(1-2.*beta_vec.^2); % why set beta to this value?
        mog_std = sqrt(1 - 2 .* beta_vec.^2 .* f_vec .* (1-f_vec));
        genotype_phenotype_prod = zeros(iters,1);
        genotype_sqr = zeros(iters,1); 
        contigency_table = zeros(iters,4); 
        for j=1:iters
            switch snp_type
                case 'binary'
                    [QTL_vec genotype_vec] = ...
                        MixtureOfGaussiansSimulateData([1-f_vec, f_vec], ...
                        [-beta_vec/2 beta_vec/2], mog_std .* [1 1], num_cases); % simulate genotype (allelic) and phenotype (QTL)
                case 'diploid'
                    [QTL_vec genotype_vec] = ...
                        MixtureOfGaussiansSimulateData([(1-f_vec)^2, 2*f_vec*(1-f_vec), f_vec^2], ...
                        [-mog_mu 0 mog_mu], mog_std .* [1 1 1], num_cases); % simulate genotype (allelic) and phenotype (QTL)
            end
            % [mog_mean mog_V] = MixtureOfGaussiansMoments([(1-f_vec)^2, 2*f_vec*(1-f_vec), f_vec^2], ...
            %    [-mog_mu 0 mog_mu], mog_std .* [1 1 1]);
            %        mog_empirical_mean = mean(QTL_vec);
            %        mog_empirical_V = var(QTL_vec);
            %        observed_beta_vec = corr(vec2column(QTL_vec), vec2column(genotype_vec)) ./ ...; % Determine empirical correlation
            %            std(genotype_vec);
            if(length(unique(genotype_vec)) == 1) % no alternative genotype seen. Effect size zero
                observed_beta_vec = 0.0000000000000000000000000000001;
            else % fit a linear regression between genotype and phenotype
                observed_beta_vec = polyfit(genotype_vec, QTL_vec, 1); observed_beta_vec = observed_beta_vec(1);
            end
            observed_beta_vec = min_abs(observed_beta_vec, 1 / sqrt(2.*f_vec.*(1-f_vec)));
            observed_f_vec = mean(genotype_vec-1) ./ 2;
            contigency_table(j,:) = QTL_params_to_p_mat(observed_f_vec, observed_beta_vec, ones(size(observed_f_vec)));        %        [observed_f_vec observed_beta_vec; 0 1]; % 2X2 format for QTL
            genotype_phenotype_prod(j) = sum(QTL_vec .* genotype_vec);
            genotype_sqr(j) = sum(genotype_vec .^ 2); % compute additional sufficient statistics
        end % loop on iters
end

