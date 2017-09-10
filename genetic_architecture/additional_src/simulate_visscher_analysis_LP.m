% Simulate parents and offspring and regress them for the LP model
%
% Input:
% n_samples - number of individuals to simulate
% N - number of pathways
% h_x - heritability of each pathway
% h_shared_env - level of shared environment (up to sibs???)
% IBD_mean - mean value of IBD (default is half for siblings)
% IBD_std - variation in IBD
% operation - just sibs (0.5) or full IBD [0,1] range
% compute_mode - numeric (compute integrals) or sampling
% trait_type - quantitative (default) or binary
% mu - prevalence (for binary trait)
% IBD_range - (optional) range at which to compute IBD (default is 0 to 1 with 0.02 resolution) 
%
% Output:
% beta - regression coefficient of response on selection strength
% IBD_sharing_vec - how much IBD is shared between individuals
% qtl_R - correlation in trait values between relatives (quantitative trait)
% lambda_R - increased in risk for relatives (bianry trait)
% power_vec - ???
% h_loci - heritability explained by all loci
% h_pop - epidemiological heritability estimate (using ACE model)
%
function [BETA IBD_sharing_vec qtl_R lambda_R ...
    power_vec h_loci h_pop] = ...
    simulate_visscher_analysis_LP(n_samples_vec, N, h_x, h_shared_env, IBD_mean, IBD_std, ...
    operation, compute_mode, trait_type, mu, IBD_range)

AssignGeneralConstants;
alpha=0.05; % significance level for computing power
iters=5000; % number of iterations for power calculations when sampling
plot_flag=0;
if(~exist('operation', 'var') || isempty(operation))
    operation = 'sibling';
end
if(~exist('compute_mode', 'var') || isempty(compute_mode))
    compute_mode = 'numeric';
end
if(~exist('trait_type', 'var') || isempty(trait_type))
    trait_type = 'quantitative';
end
if(~exist('h_shared_env', 'var') || isempty(h_shared_env))
    h_shared_env = 0;
end

n_samples = max(n_samples_vec); % compute using the maximum sample size
switch trait_type
    case {'binary', 'disease'}
        options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance
        K=1; mu_l = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu)), 0, 1, options); % find mu_l that keeps the prevalence
        x_mu_l = norminv(1-mu_l); % set threshold for disease
end

switch operation
    case 'sibling' % compute just variation in IBD for siblings
        if(~exist('IBD_mean', 'var') || isempty(IBD_mean))
            IBD_mean = 0.5;
        end
        if(~exist('IBD_std', 'var') || isempty(IBD_std))
            IBD_std = IBD_mean/10;
        end
        switch compute_mode
            case 'numeric'
                num_IBD=3; IBD_vec = [0.49 0.5 0.51];
            case 'sampling'
                num_IBD=1;
        end
    case 'full-range' % compute the entire range of IBD from 0 to 1 
        if(~exist('IBD_range', 'var') || isempty(IBD_range))
            res = 0.02;
            IBD_range = linspace(0, 1, 1+1/res);
        end
        % New: compute average phenotypic correlation for ALL IBD
        num_IBD=length(IBD_range);
end
use_IBD_std = IBD_std * sqrt(N); % adjust st.d. such that overall genotype has correct IBD variance


%qtl_R = []; lambda_R = [];
switch trait_type % compute h_loci, h_pop and relatives correlation
    case {'quantitative', 'QTL'}
        [~, ~, ~, ~, ~, h_loci, tmp_h_pop] = ...
            compute_k_of_N_gaussian_statistics([0 0], [1 1], h_x, h_shared_env, [], ...
            'MAX', [], N, 10, 'numeric', 0, {'ACE'});
        h_pop  = tmp_h_pop{1};
    case {'binary', 'disease'} % need to compute h_loci and h_pop from here!!!
        [lambda_R STATS] = ...             % Need to extract h_loci from here!
            compute_k_of_N_liabilities_statistics(N, 1, mu_l, h_x, h_shared_env, ...
            2);
        h_loci = STATS.h_liab_loci;
        h_pop = STATS.h_liab_twins;
end % switch

IBD_sharing_vec = zeros(num_IBD,1);
qtl_R = zeros(num_IBD,1);  lambda_R = qtl_R;
for IBD_ind = 1:num_IBD % compute familial correlation/risk
    switch operation
        case 'full-range'
            sprintf('run_IBD %ld out of %ld', IBD_ind, num_IBD)
            IBD_vec = repmat(IBD_range(IBD_ind), n_samples, N);
            IBD_sharing_vec(IBD_ind) = IBD_range(IBD_ind);            
            switch trait_type
                case {'quantitative', 'QTL'}
                    [~, tmp_qtl_R] = ...
                        qtl_familial_correlations_internal(N, N, IBD_range(IBD_ind), ...
                        'MAX', [], h_x, h_shared_env, 'numeric', []);
                    qtl_R(IBD_ind) = tmp_qtl_R(1);
                case {'binary', 'disease'}
                    tmp_lambda_R = compute_lambda_R_LP_internal(N, 1, ...
                        h_x, h_shared_env, IBD_range(IBD_ind), mu_l, mu);
                    lambda_R(IBD_ind) = tmp_lambda_R(1);
            end
        case 'sibling'
            switch compute_mode
                case 'numeric'
                    IBD_sharing_vec(IBD_ind) = IBD_vec(IBD_ind);
                    switch trait_type
                        case {'quantitative', 'QTL'}
                            [~, tmp_qtl_R] = ...
                                qtl_familial_correlations_internal(N, N, IBD_vec(IBD_ind), ...
                                'MAX', [], h_x, h_shared_env, 'numeric', []);
                            qtl_R(IBD_ind) = tmp_qtl_R(1);
                            %                            [mu_z sigma_z] = maxnormstat(N);
                        case {'binary', 'disease'}
                            tmp_lambda_R = compute_lambda_R_LP_internal(N, 1, ...
                                h_x, h_shared_env, IBD_vec(IBD_ind), mu_l, mu);
                            lambda_R(IBD_ind) = tmp_lambda_R(1);
                    end % switch trait type
                    
                case 'sampling' % sample
                    power_mat = zeros(length(n_samples_vec),iters);
                    for cur_iter = 1:iters % loop to compute power empirically
                        if(mod(cur_iter, 50) == 0)
                            run_iter = cur_iter
                        end
                        IBD_vec = randn(n_samples,N) .* use_IBD_std + IBD_mean; % IBD for each pathway for each pair of siblings
                        IBD_vec = min(max(IBD_vec,0), 1);
                        %    end
                        shared_genotype_vec = randn(n_samples,N).*sqrt(h_x) .* sqrt(IBD_vec); % sqrt(IBD_vec);
                        for i=1:2 % loop on two siblings
                            unique_genotype_vec{i} = randn(n_samples,N).* ...
                                sqrt(h_x) .* sqrt(1-IBD_vec); % sqrt(1-IBD_vec);
                            unique_environmental_vec{i} = randn(n_samples,N).*sqrt(1-h_x);
                            pathway_vec{i} = unique_genotype_vec{i} + ...
                                unique_environmental_vec{i}+shared_genotype_vec;
                            phenotype_qtl_vec(:,i) = max(pathway_vec{i},[],2);
                            switch trait_type
                                case {'binary', 'disease'} % take only if it's above a threshold
                                    phenotype_qtl_vec(:,i) = (phenotype_qtl_vec(:,i)>x_mu_l);
                            end
                            
                        end
                        mu_z = mean(phenotype_qtl_vec(:));
                        sigma_z = 0.5*(std(phenotype_qtl_vec(:,1)) + std(phenotype_qtl_vec(:,2)));
                        
                        % corr_phenotype_R = (phenotype_qtl_vec(:,1)-phenotype_qtl_vec(:,2)).^2;
                        qtl_R = ... % new: compute correlation coefficient!!!
                            (phenotype_qtl_vec(:,1)-mu_z).*(phenotype_qtl_vec(:,2)-mu_z) ./ ...
                            sigma_z^2;
                        IBD_sharing_vec = sum(IBD_vec,2) ./ N;
                        for n=1:length(n_samples_vec)
                            [TEMP_BETA BETA_INTERVAL] = ...
                                regress(qtl_R(1:n_samples_vec(n)), ...
                                [IBD_sharing_vec(1:n_samples_vec(n)) ...
                                ones(n_samples_vec(n), 1)], alpha);
                            %BETA_MAT(cur_iter,n) = TEMP_BETA(1);
                            if((h_pop <  BETA_INTERVAL(1,1)) || (h_pop > BETA_INTERVAL(1,2)))
                                power_mat(n,cur_iter) = 1;
                            end
                        end % loop on sample size
                        
                        %                     mean_corr_phenotype_R(IBD_ind) = mean(corr_phenotype_R);
                        %                     mean_IBD_sharing_vec(IBD_ind) = mean(IBD_sharing_vec);
                        %                     corr_corr_phenotype_R(IBD_ind) = corr(phenotype_qtl_vec(:,1), phenotype_qtl_vec(:,2));
                        
                    end % loop on iterations for power computation
                    power_vec = mean(power_mat, 2);
            end % switch compute mode (sampling/numeric integral)
    end % switch operation (sibling or full range IBD)
    switch trait_type
        case 'quantitative'
            
        case 'binary' % transform from lambda_R to r_R (qtl_R)
            qtl_R(IBD_ind) = ...
                familial_risk_to_heritability(lambda_R(IBD_ind), 'liability', mu, 1);
%            lambda_R(IBD_ind) = qtl_R(IBD_ind);
    end    
end % loop on different IBD_vecs


%figure;
if(plot_flag)
    switch operation
        case 'sibling'
            cur_color='b'; cur_symbol = '.';
        case 'full-range'
            cur_color=color_vec(N); cur_symbol = '-';
            %        corr_phenotype_R  = corr_corr_phenotype_R; % mean_ ? what to plot for phenotype?
            %        IBD_sharing_vec = mean_IBD_sharing_vec;
    end
    hold on; plot(IBD_sharing_vec, qtl_R, [cur_symbol cur_color]);
    %    rho = corr(vec2column(IBD_sharing_vec), vec2column(qtl_R))
end


switch operation
    case 'sibling'
        [BETA] = polyfit(vec2column(IBD_sharing_vec), vec2column(qtl_R), 1);
        [BETA_AGAIN BETA_INTERVAL] = regress(vec2column(qtl_R), ...
            [vec2column(IBD_sharing_vec) ones(length(IBD_sharing_vec), 1)], alpha);
        if(plot_flag)
            plot(IBD_sharing_vec, IBD_sharing_vec.*BETA(1)+BETA(2), 'r');
            title(['Regression of \Delta genotype vs. \Delta phenotype LP(N=' ...
                num2str(N) ', h_x=' num2str(h_x*100,2) '%)' ...
                '% \beta=' num2str(BETA(1)*100,3) '% -\beta/2=' num2str(-BETA(1)*100/2,2) '% h_{loci}=' ... % ' \rho=' num2str(rho*100,2)
                num2str(h_loci*100,2) '% r_{DZ}=' num2str(qtl_R(2)*100,2) '%']);
        end
        
        switch compute_mode       % Compute power here
            case 'numeric'
                x_alpha = norminv(alpha); % perform one-sided LEFT test
                T_stat = BETA(1) - min(1,h_pop); % this is the noncentrality parameter
                T_stat = T_stat * IBD_std /1;
                power_vec = normcdf(x_alpha-T_stat.*sqrt(n_samples_vec));
        end
        BETA = BETA(1); 
    case 'full-range'
        mid_ind = ceil(0.5*num_IBD); % size(qtl_R, 2));
        %        BETA = polyfit(IBD_sharing_vec([mid_ind-1 mid_ind+1]), ...
        %            qtl_R([mid_ind-1 mid_ind+1]), 1);
        BETA = (qtl_R(mid_ind-1)- qtl_R(mid_ind+1)) / ...
            (IBD_sharing_vec(mid_ind-1)- IBD_sharing_vec(mid_ind+1));
        
        if(plot_flag)
            title(['Regression of \Delta genotype vs. \Delta phenotype LP(N=' ...
                num2str(N) ', h_x=' num2str(h_x*100,2) '%)' ...
                '% h_{loci}=' ... % ' \rho=' num2str(rho*100,2)
                num2str(h_loci*100,2) '% r_{DZ}=' num2str(qtl_R(2)*100,2) '%']);
            %        plot(IBD_sharing_vec, corr_qtl_R_analytic, cur_color); % plot analytic
        end
end % switch operation

if(plot_flag)
    xlabel('IBD (%)');
    ylabel('corr(Z_1,Z_2)'); % ylabel('(\Delta Z)^2');
end

IBD_sharing_vec = vec2row(IBD_sharing_vec);
qtl_R = vec2row(qtl_R);

if(~exist('power_vec', 'var'))
    power_vec = [];
end


