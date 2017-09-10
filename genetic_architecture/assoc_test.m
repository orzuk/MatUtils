% Perform test for association
%
% Input:
% contigency_table_vec - data counts (1x4, 1x6 or 1x8 table) of genotype/allele and trait
% test_type - type of signal to test (single locus, epistasis etc.)
% test_stat - statistic to used for testing
% n_samples_vec - number of samples (case? controls?)
% full_flag - flag saying if simulate full set of genotype/phenotype or just summary statistics
% plot_stat_figures - flag saying if to plot figures
% model_params - parameters of null model
%
% Output:
% p_vals_vec - p-value of test
% stat_vec - statistic used (usually chi-square, but can be also t-test statistic)
%
function [p_vals_vec stat_vec] = assoc_test(...
    contigency_table_vec, test_type, test_stat, n_samples_vec, ...
    full_flag, plot_stat_figures, model_params, alpha_vec) % Perform association test

s = warning('off', 'all');

if(~exist('full_flag', 'var') || isempty(full_flag))
    full_flag = 0; % set full_flag (default: NO)
end
if(~exist('n_samples_vec', 'var') || isempty(n_samples_vec))
    n_samples_vec = sum(contigency_table_vec,2); % set sample sizes. Works only for binary traits
end
num_cells = size(contigency_table_vec, 2); % number of different possiblies (can be 4,6,8)
iters = size(contigency_table_vec, 1); % different number of effects
if(length(n_samples_vec) == 1)
    n_samples_vec = repmat(n_samples_vec, iters, 1);
end
%n_cases_vec = sum(contigency_table_vec(:,1:num_cells/2),2); % they're used only in armitage - where they're specified
%n_controls_vec = sum(contigency_table_vec(:,(num_cells/2+1):num_cells),2);


switch test_stat % perform test
    case 'hypergeometric' % assume 2x2 tables
        p_vals_vec = 1-hygecdf(max(0,contigency_table_vec(:,4)-1), n_samples_vec, ...
            sum(contigency_table_vec(:,[2 4]),2), ...
            sum(contigency_table_vec(:,[3 4]),2));
    case 'chi-square'
        switch lower(test_type)
            case {'single-locus', 'marginal'} % this test treats SNPs as binary
                observed_table_vec = contigency_table_vec;
                p1_vec = sum(contigency_table_vec(:,[2 4]),2) ./ n_samples_vec; % Pr(x_1 = 1)
                p2_vec = sum(contigency_table_vec(:,[3 4]),2) ./ n_samples_vec; % Pr(x_2 = 1) or Pr(Z = 1)
                expected_table_vec = [(1-p1_vec) .* (1-p2_vec) ...
                    p1_vec .* (1-p2_vec)  (1-p1_vec).*p2_vec p1_vec.*p2_vec] .* ...
                    repmat(n_samples_vec, 1, 4);
            case {'armitage', 'additive', 'trend'} % trend test (for 2x3 test)
                table_format = 'interleaved';
                switch table_format % how does the input contigency table look
                    case 'controls_cases'
                        n_controls_vec = sum(contigency_table_vec(:,1:num_cells/2),2); % first controls appearing, then cases !
                        n_cases_vec = sum(contigency_table_vec(:,(num_cells/2+1):num_cells),2); % flipped convension
                        
                        observed_table_vec = [2.*contigency_table_vec(:,1) + contigency_table_vec(:,2), ...
                            2.*contigency_table_vec(:,3) + contigency_table_vec(:,2), ...
                            2.*contigency_table_vec(:,4) + contigency_table_vec(:,5), ...
                            2.*contigency_table_vec(:,6) + contigency_table_vec(:,5)] ./ 2; % collapse to allelic table
                        p1_vec = ( sum(contigency_table_vec(:,[2 5]),2) + ... % one allele
                            2.*sum(contigency_table_vec(:,[3 6]),2) ) ./ (2*n_samples_vec); % Pr(x_1 = 1)
                        p2_vec = sum(contigency_table_vec(:,[4:6]),2) ./ n_samples_vec; % Pr(Z = 1)
                        expected_table_vec = [(1-p1_vec) .* (1-p2_vec) ...
                            p1_vec .* (1-p2_vec)  (1-p1_vec).*p2_vec p1_vec.*p2_vec] .* ...
                            repmat(n_samples_vec, 1, 4); % expected stays the same (assumes HW equilibrium)
                        
                        Y_sqr = n_samples_vec .* (n_samples_vec.*(2.*contigency_table_vec(:,6) + contigency_table_vec(:,5)) - ...
                            n_cases_vec.*(2.*sum(contigency_table_vec(:,[3 6]),2) + sum(contigency_table_vec(:,[2 5]),2))).^2;   % compute Armitage test also directly
                        Y_sqr_denom = n_cases_vec.*n_controls_vec .* (n_samples_vec.*(4.*sum(contigency_table_vec(:,[3 6]),2) + ...
                            sum(contigency_table_vec(:,[2 5]),2)) - ...
                            (2.*sum(contigency_table_vec(:,[3 6]),2) + sum(contigency_table_vec(:,[2 5]),2)).^2);
                    case 'interleaved' % different format
                        n_controls_vec = sum(contigency_table_vec(:,1:2:end),2); % first controls appearing, then cases !
                        n_cases_vec = sum(contigency_table_vec(:,2:2:end),2); % flipped convension
                        
                        observed_table_vec = [2.*contigency_table_vec(:,1) + contigency_table_vec(:,3), ...
                            2.*contigency_table_vec(:,5) + contigency_table_vec(:,3), ...
                            2.*contigency_table_vec(:,2) + contigency_table_vec(:,4), ...
                            2.*contigency_table_vec(:,6) + contigency_table_vec(:,4)] ./ 2; % collapse to allelic table
                        p1_vec = ( sum(contigency_table_vec(:,[3 4]),2) + ... % one allele
                            2.*sum(contigency_table_vec(:,[5 6]),2) ) ./ (2*n_samples_vec); % Pr(x_1 = 1)
                        p2_vec = sum(contigency_table_vec(:,[2 4 6]),2) ./ n_samples_vec; % Pr(Z = 1)
                        expected_table_vec = [(1-p1_vec) .* (1-p2_vec) ...
                            p1_vec .* (1-p2_vec)  (1-p1_vec).*p2_vec p1_vec.*p2_vec] .* ...
                            repmat(n_samples_vec, 1, 4); % expected stays the same (assumes HW equilibrium)
                        
                        Y_sqr = n_samples_vec .* (n_samples_vec.*(2.*contigency_table_vec(:,6) + contigency_table_vec(:,4)) - ...
                            n_cases_vec.*(2.*sum(contigency_table_vec(:,[5 6]),2) + sum(contigency_table_vec(:,[3 4]),2))).^2;   % compute Armitage test also directly
                        Y_sqr_denom = n_cases_vec.*n_controls_vec .* (n_samples_vec.*(4.*sum(contigency_table_vec(:,[5 6]),2) + ...
                            sum(contigency_table_vec(:,[3 4]),2)) - ...
                            (2.*sum(contigency_table_vec(:,[5 6]),2) + sum(contigency_table_vec(:,[3 4]),2)).^2);
                end
                
                Y_sqr = Y_sqr ./ Y_sqr_denom;
                
                
            case {'pairwise', 'epistasis'} % compute chi-square for the interaction to be significant
                % We assume that the null model is: Pr(Z=1|x_i,x_j) =
                % \alpha_1 * x_1 + \alpha_2 * x_2 + \beta.
                % The test sums Pr(00*)+Pr(11*) and compares them to
                % Pr(01*)+Pr(11*) and then does a standard 2x2
                % chi-square with one deg. freedom
                % %                     empirical_p_x_12_vec = contigency_table_vec(:,1:2:7) + ...
                % %                         contigency_table_vec(:,2:2:8); % Pr(x_1,x_2) all four possibilites
                % %                     empirical_p_x_1_vec = [sum(contigency_table_vec(:,1:4)) ...
                % %                         sum(contigency_table_vec(:,5:8))];
                % %                     empirical_p_x_2_vec = [sum(contigency_table_vec(:,[1 2 5 6])) ...
                % %                         sum(contigency_table_vec(:,[3 4 7 8]))];
                % %                     empirical_p_z_given_x_12_vec = contigency_table_vec(:,2:2:8) ./ ...
                % %                         (contigency_table_vec(:,1:2:7) + contigency_table_vec(:,2:2:8)); % Pr(z=1|x_1,x_2)
                % %                     empirical_p_z_given_x_1_vec = (contigency_table_vec(:,2:2:8) + ...
                % %                         contigency_table_vec(:,2:2:8)) ./ ...
                % %                         (contigency_table_vec(:,2:2:8) + contigency_table_vec(:,2:2:8));
                
                empirical_p_z_vec = sum(contigency_table_vec(:,2:2:8),2) ./ ...
                    n_samples_vec; % Pr(z = 1)
                empirical_p_x_1_xor_x2_vec = 1-sum(contigency_table_vec(:,[1 2 7 8]),2) ./ ...
                    n_samples_vec; % Pr(x_1 + x_2 = 0 mod 2)
                expected_table_vec = ... % Pr(z) * Pr(x_1+x_2)
                    [(1-empirical_p_z_vec) .* (1-empirical_p_x_1_xor_x2_vec) ...
                    empirical_p_z_vec .* (1-empirical_p_x_1_xor_x2_vec) ...
                    (1-empirical_p_z_vec).*empirical_p_x_1_xor_x2_vec ...
                    empirical_p_z_vec.*empirical_p_x_1_xor_x2_vec] .* n_samples_vec;
                
                %                     expected_table_vec = ... % Pr(z) * Pr(x_1+x_2)
                %                         0.5 .* [(1-empirical_p_x_1_xor_x2_vec) ...
                %                          (1-empirical_p_x_1_xor_x2_vec) ...
                %                         empirical_p_x_1_xor_x2_vec ...
                %                         empirical_p_x_1_xor_x2_vec] .* n_samples_vec(i);
                observed_table_vec = ...  % P(00*)+P(11*), P(01*)+P(10*)
                    [contigency_table_vec(:,1) + contigency_table_vec(:,7) ...
                    contigency_table_vec(:,2) + contigency_table_vec(:,8) ...
                    contigency_table_vec(:,3) + contigency_table_vec(:,5) ...
                    contigency_table_vec(:,4) + contigency_table_vec(:,6)];
                
                % % % %                 theoretical_table_vec = [p_vec(1) + p_vec(7) ...
                % % % %                     p_vec(2) + p_vec(8) ...
                % % % %                     p_vec(3) + p_vec(5) ...
                % % % %                     p_vec(4) + p_vec(6)] .* n_samples_vec;
            case {'likelihood-ratio-test', 'llr'} % problem: computing MLE might be hard for certain no-epistasis models!!!!
                
        end % switch test type
        chi_stat_vec = sum( (observed_table_vec - expected_table_vec).^2 ./ ...
            expected_table_vec, 2 );
        if(~exist('Y_sqr', 'var'))
            stat_vec = chi_stat_vec; % copy to output
        else
            stat_vec = Y_sqr;
        end
        p_vals_vec = 1-chi2cdf(stat_vec, 1); % one degree of freedom
        
    case {'probit', 'logit', 'logistic', 'probit-chi-square'} % just fit a logistic regression and get p-value of each term.
        % Fitting formula is: Pr(z = 1 | x_1, x_2) / Pr(z = 0 | x_1, x_2) =
        % log(alpha_vec + alpha_1*x_1 + alpha_2*x_2 + alpha_12*x_1*x_2.
        % This is the 'full' (ravui) model - should be fast to fit!!!
        adjust_case_control_to_population=0; % not working yet
        if(full_flag)
            %            x_12_logit_matrix = [contigency_table_vec(:,:,1) contigency_table_vec(:,:,2) ...
            %                contigency_table_vec(:,:,1) .* contigency_table_vec(:,:,2)];
            z_logit_vec = contigency_table_vec(:,:,3)';
            N_logit_vec = [];
            optional_inds = vec2row(find(min(z_logit_vec) < max(z_logit_vec))); % require at least one count for each (x_i,x_j)
        else % work on summary statistics (faster, more complicated)            
            x_12_logit_matrix = [0 0 0; 0 1 0; 1 0 0; 1 1 1]; % x_1, x_2, x_1&x_2
            if(adjust_case_control_to_population)
                mu = model_params(3);
                contigency_table_vec_population = ...
                    case_control_prob_to_pop_prob( ...
                    contigency_table_vec./n_samples_vec(1), mu).*n_samples_vec(1);
                z_logit_vec = contigency_table_vec_population(:, (2:2:8))'; % counts of Z equal one
                N_logit_vec = (contigency_table_vec_population(:, (1:2:7)) ...
                    + contigency_table_vec_population(:, (2:2:8)))'; % counts of Z equal one and zero
            else
                z_logit_vec = contigency_table_vec(:, (2:2:8))'; % counts of Z equal one
                N_logit_vec = (contigency_table_vec(:, (1:2:7)) ...
                    + contigency_table_vec(:, (2:2:8)))'; % counts of Z equal one and zero
            end
            
            optional_inds = vec2row(find(min(z_logit_vec))); % require at least one count for each (x_i,x_j)
            %%            z_logit_vec = [z_logit_vec' N_logit_vec']';
        end
        stat_vec = zeros(iters,1); p_vals_vec = ones(iters,1); % just put ones (no power at all) as default
        for j=optional_inds % must loop ... this part is slow
            if(mod(j,50)==0)
                sprintf('n=%ld, run_iter %ld ', n_samples_vec(1), j)
            end
            if(full_flag) % here cur_z_vec holds ..
                x_12_logit_matrix = [contigency_table_vec(j,:,1)' contigency_table_vec(j,:,2)' ...
                    contigency_table_vec(j,:,1)' .* contigency_table_vec(j,:,2)'];
                cur_z_vec = z_logit_vec(:,j);
            else % here we also hold N
                cur_z_vec = [z_logit_vec(:,j) N_logit_vec(:,j)];
            end
            %            warning('off', 'MATLAB:glmfit'); % for some reason this doesn't work
            [~, ~, stats] =  glmfit(x_12_logit_matrix, ... % allow both logistic and probit regression
                cur_z_vec, 'binomial', 'link', ...
                strrep( str2word('-', test_stat, 1), 'stic', 't'));   % [z_logit_vec(:,j) N_logit_vec(:,j)],
            z_ravui_probs_vec = normcdf(stats.beta(1) + ...
                stats.beta(2) .* x_12_logit_matrix(:,1) + ...
                stats.beta(3) .* x_12_logit_matrix(:,2) + ...
                stats.beta(4) .* x_12_logit_matrix(:,3));
            observed_z_probs_vec = cur_z_vec(:,1) ./ cur_z_vec(:,2);
            
            
            switch test_stat
                case 'probit-chi-square' % special case: perform chi-square to detect epistasis
                    [~, ~, stats_additive] =  glmfit(x_12_logit_matrix(:,1:2), ... % allow both logistic and probit regression
                        cur_z_vec, 'binomial', 'link', ...
                        strrep( str2word('-', test_stat, 1), 'stic', 't'));   % [z_logit_vec(:,j) N_logit_vec(:,j)],
                    z_probs_vec = normcdf(stats_additive.beta(1) + ...
                        stats_additive.beta(2) .* x_12_logit_matrix(:,1) + ...
                        stats_additive.beta(3) .* x_12_logit_matrix(:,2));
                    input_h_x = model_params(1);
                    freq = model_params(5);
                    theoretical_beta(2) = sqrt(input_h_x / (freq .*(1-freq)));
                    theoretical_beta(3) = theoretical_beta(2);
                    theoretical_beta(1) = -2*theoretical_beta(2) * freq;
                    theoretical_beta = theoretical_beta ./ sqrt(1-2*input_h_x);
                    if(~isreal(theoretical_beta))
                        XXXXX = 1234124123;
                    end
                    z_theoretical_probs_vec = normcdf(theoretical_beta(1) + ...
                        theoretical_beta(2) .* x_12_logit_matrix(:,1) + ...
                        theoretical_beta(3) .* x_12_logit_matrix(:,2));
                    
                    %                    observed_z_probs_vec = cur_z_vec(:,1) ./ cur_z_vec(:,2)
                    tmp_expected_table = contigency_table_vec(j,1:2:end) + contigency_table_vec(j,2:2:end); % compute expected table
                    expected_table(j,:) =  ...
                        mat2vec([tmp_expected_table' .* (1-z_probs_vec) ...
                        tmp_expected_table' .* (z_probs_vec)]')';
                    if(adjust_case_control_to_population)
                        expected_table(j,:) = ...
                            pop_prob_to_case_control_prob(expected_table(j,:));
                    end
                    observed_table(j,:) = contigency_table_vec(j,:);
                    stat_vec(j) = sum((expected_table(j,:) - observed_table(j,:)).^2 ./ ...
                        expected_table(j,:)); % compute chi-square statistics with 1 d.f.
                    p_vals_vec(j) = 1-chi2cdf(stat_vec(j),1);
                otherwise
                    stat_vec(j) = stats.t(4); % t statistic for interaction term
                    p_vals_vec(j) = stats.p(4); % p-value of interaction term
            end
        end % loop on samples which were not constant
        
        
    case {'goodness-of-fit', 'LRT', 'likelihood-ratio-test', ...
            'gof2d', 'goodness-of-fit-two-dim', 'gof-squares'} % Perform a chi-square(?) test. Test expected vs. observed
        if(exist('model_params', 'var') && (~isempty(model_params)))
            input_h_x = model_params(1); input_h_x_MLT = model_params(2);
            fit_h_x = 0; % Assume we get h_x as input as part of the model
        else
            fit_h_x = 1;
        end
        stat_vec = zeros(iters,1); p_vals_vec = ones(iters,1);
        if(full_flag)
            % N=1; % We assume that we've got ONE liability and both pathways are in it
            
            h_x = zeros(iters, 1);
            optional_inds = vec2row(find(min(contigency_table_vec(:,:,3),[],2) < ...
                max(contigency_table_vec(:,:,3),[],2))); % require at least one count for each (x_i,x_j)
            
            if(fit_h_x) % we can fit this parameter or get as input
                mu = mean(mat2vec(contigency_table_vec(:,:,3))); % prevalence of DISEASE !!!
                for j=optional_inds % 1:iters
                    h_x(j) = corr(contigency_table_vec(j,:,1)', contigency_table_vec(j,:,3)')^2 + ...
                        corr(contigency_table_vec(j,:,2)', contigency_table_vec(j,:,3)')^2;
                end
                % %              for j=optional_inds % 1:iters
                % %                 h_x2(j) = corr(contigency_table_vec(j,:,1)'+contigency_table_vec(j,:,2)', ...
                % %                     contigency_table_vec(j,:,3)')^2;
                % % %                + ...
                % % %                    corr(contigency_table_vec(j,:,2)', contigency_table_vec(j,:,3)')^2;
                % %             end
                h_x = min(heritability_scale_change(h_x, 'liability', mu), 1); % This is SINGLE-LT heritability
            else
                mu = model_params(3); % input also the prevalence
                h_x(:) = 2*input_h_x; % account for two X's
            end
            x_mu = norminv(1-mu); % threshold assuming single LT model
            
            if(~exist('plot_stat_figures', 'var') || isempty(plot_stat_figures))
                plot_stat_figures = 0;
            end
            
            if( (plot_stat_figures) || ismember(test_stat, {'LRT', 'likelihood-ratio-test'}) )% Plot statistic and see tat it matches chi square !!!
                N=2;
                mu_l = fminbnd(@(x) abs(binocdf(1-1, 2, x)-(1-mu)), 0, 1); % find mu_l that keeps the PREVALENCE
                x_mu_l = norminv(1-mu_l);
                h_x_MLT = h_x;
                if(fit_h_x)
                    for j=1:iters
                        convert_j_h_x_to_MLT = j
                        h_x_MLT(j) = heritability_scale_change_MLT(h_x(j), 1, 2, mu, 'MLT'); % compute assiuming (2,1) model. Take a long time
                    end
                else
                    h_x_MLT(:) = input_h_x_MLT; % heritability_scale_change_MLT(h_x(1), 1, 2, mu, 'MLT'); % compute assiuming (2,1) model. Take a long time
                end
                z_expected_MLT = 1 - ...
                    (1-mu_l)^(N-2) .* normcdf( (x_mu_l - contigency_table_vec(:,:,1) .* ...
                    repmat(sqrt(h_x_MLT), 1, n_samples_vec(1))) ./ ...
                    sqrt(1-repmat(h_x_MLT, 1, n_samples_vec(1))) ) .* ...
                    normcdf( (x_mu_l - contigency_table_vec(:,:,2) .* ...
                    repmat(sqrt(h_x_MLT), 1, n_samples_vec(1))) ./ ...
                    sqrt(1-repmat(h_x_MLT, 1, n_samples_vec(1))) ); % This is the expected value under MLT model
            end % if plot figures
            % % %             for j=1:2 % Normalize
            % % %                 contigency_table_vec(:,:,j) = repmat(sqrt(h_x./2), 1, n_samples_vec(1)) .* ...
            % % %                     contigency_table_vec(:,:,j) ./ ...
            % % %                     repmat(std(contigency_table_vec(:,:,1),[],2), 1, n_samples_vec(1));
            % % %             end
            z_expected =  1 - ...  % Expectation of Z under the null model (LT)
                normcdf( (x_mu - (contigency_table_vec(:,:,1) + contigency_table_vec(:,:,2)) .* ...
                repmat(sqrt(0.5*h_x), 1, n_samples_vec(1))  ) ./ ...
                sqrt(1-repmat(h_x, 1, n_samples_vec(1)))); % h_x is already the sum of the heritability for two variables
            
            %.* ...
            %    normcdf( (x_mu - contigency_table_vec(:,:,2)) ./ sqrt(1-h_x)); % Expected value under liability treshold model
            
            switch test_stat
                case {'goodness-of-fit'}
                    z_diff = contigency_table_vec(:,:,3) - z_expected;
                    stat_vec = sum(z_diff.^2 ./ z_expected, 2);
                    stat_vec(setdiff(1:iters, optional_inds)) = 0; % don't
                    
                    stat_mu = n_samples_vec(1) - sum(z_expected,2); % Compute statistic moment
                    stat_sigma = sqrt(  sum((1-2.*z_expected).^2.*(1-z_expected) ./ z_expected,2) );
                    stat_vec = (stat_vec - stat_mu) ./ stat_sigma;
                    p_vals_vec = 1-normcdf(stat_vec); % Give p-values, based on the same
                    test_params = [];
                case {'LRT', 'likelihood-ratio-test'} % Don't normalize !!!
                    stat_vec = sum(log (1 - contigency_table_vec(:,:,3) + ...
                        z_expected_MLT .* (2.*contigency_table_vec(:,:,3)-1)),2) - ...
                        sum(log (1 - contigency_table_vec(:,:,3) + ...
                        z_expected .* (2.*contigency_table_vec(:,:,3)-1)),2);
                    SSS = log (1 - contigency_table_vec(:,:,3) + ...
                        z_expected_MLT .* (2.*contigency_table_vec(:,:,3)-1)) - ...
                        log (1 - contigency_table_vec(:,:,3) + ...
                        z_expected .* (2.*contigency_table_vec(:,:,3)-1));
                    
                    
                    p_vals_vec = 1 ./ (1+exp(stat_vec)); % What's the LRT distribution? (no parameters fitted here!). Use a bayesian prior with (0.5,0.5) s
                    test_params = [];
                    [LRT_mu LRT_var] = ...
                        LRT_stat_moments_MLT(2, 1, mu, h_x(1)/2, 'MLT'); % here h_x is on the liability scale
                    test_params = [LRT_mu LRT_var];
                    
                case {'gof2d', 'goodness-of-fit-two-dim', 'gof-squares'} % divide the (X1,X2) space to equiprobable bins
                    points_in_bin = 50; % This is the average number of points per bin
                    num_bins = max(2, floor(sqrt(n_samples_vec(1)/points_in_bin))); % number of bins in each axis
                    x_grid = [-100 norminv((1:num_bins-1) ./ num_bins) 100];
                    z_expected_table = zeros(num_bins,num_bins,iters);
                    z_observed_table = zeros(num_bins,num_bins,iters);
                    z_counts_table = zeros(num_bins,num_bins,iters);
                    for i=1:num_bins % loop on bins for X1
                        for j=1:num_bins % loop on bins for X2
                            good_inds_mat = zeros(iters, n_samples_vec(1));
                            good_inds = intersect(find((contigency_table_vec(:,:,1) <= x_grid(i+1)) & ...
                                (contigency_table_vec(:,:,1) >= x_grid(i))), ...
                                find((contigency_table_vec(:,:,2) <= x_grid(j+1)) & ...
                                (contigency_table_vec(:,:,2) >= x_grid(j))));
                            good_inds_mat(good_inds) = 1;
                            z_expected_table(i,j,:) = sum(z_expected .* good_inds_mat,2);
                            z_observed_table(i,j,:) = sum(contigency_table_vec(:,:,3) .* ...
                                good_inds_mat,2);
                            z_counts_table(i,j,:) = sum(good_inds_mat,2); % how many counts in each cell
                        end
                    end
                    stat_vec = reshape(sum(sum(z_counts_table.*(z_expected_table-z_observed_table).^2 ./ ...
                        (z_expected_table .* (z_counts_table-z_expected_table)))), iters, 1);
                    dof = num_bins.^2; % set numbers of degrees of freedom
                    p_vals_vec = 1-chi2cdf(stat_vec, dof);
                    
                    
                    if(plot_stat_figures) % Plot statistic and see that it matches chi square !!!% Compute the theoretical distribution (temp, remove!!!!)
                        z_expected_table_LT = zeros(num_bins,num_bins);
                        z_expected_table_MLT = zeros(num_bins,num_bins);
                        N=2;       mu = model_params(3); h_x = model_params(1:2);
                        mu_l = fminbnd(@(x) abs(binocdf(1-1, N, x)-(1-mu)), 0, 1); % find mu_l that keeps the PREVALENCE
                        x_mu = norminv(1-mu); x_mu_l = norminv(1-mu_l);
                        x_grid = [-100 norminv((1:num_bins-1) ./ num_bins) 100];
                        num_points_in_bin = round(n_samples_vec(1) / num_bins^2); % compute actual number of points
                        
                        for k=1:num_bins % loop on bins for X1
                            for j=1:num_bins % loop on bins for X2
                                z_expected_table_LT(k,j) = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
                                    z_expected_given_two_x_MLT(x1, x2, 1, h_x(1), mu, x_mu), ...% Compute expected variance
                                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)) / ...
                                    quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2), ...
                                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)); % normalize by density
                                z_expected_table_MLT(k,j) = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
                                    z_expected_given_two_x_MLT(x1, x2, N, h_x(2), mu_l, x_mu_l), ...% Compute expected variance
                                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)) / ...
                                    quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2), ...
                                    x_grid(k), x_grid(k+1), x_grid(j), x_grid(j+1)); % normalize by density
                            end
                        end
                        non_centrality_parameter = sum(sum( ...
                            num_points_in_bin.*(z_expected_table_MLT - z_expected_table_LT).^2 ./ ...
                            (z_expected_table_LT.*(1 - z_expected_table_LT)) )); % here expected are probabilities (not counts)
                        test_params(1) = dof;
                        test_params(2) = non_centrality_parameter;
                    end
                    
            end
            
            if(plot_stat_figures) % Plot statistic and see that it matches chi square !!!
                plot_stat_figures_internal(stat_vec, z_expected, z_expected_MLT, ...
                    contigency_table_vec, n_samples_vec, mu, h_x, test_stat, test_params);
            end % plot_stat_figures
            
        else % don't simulate fully
            
        end % if full flag
        
    case 'z-score' % What here?
        
    case {'QTL', 'chi-square-QTL'} % Perform an empirical test of QTL
        [f_vec beta V] = p_mat_to_QTL_params(contigency_table_vec); % encoding convension: last value is traits variance
        %        genotype_phenotype_prod =
        %        genotype_sqr = (2.*f_vec.^2 - 2.*f_vec + 1) .* n_samples_vec;
        % encoding convension:
        
        %        V_explained = beta_to_variance_explained(beta, f_vec, V, 'diploid'); % multiply effect by two
        stat_vec = zeros(iters,1); p_vals_vec = zeros(iters,1);
        for j=1:iters% Simulate a population and their value - WHY simulate? shouldn't this be as input?
            %            [cur_samples_QTL_vec cur_samples_genotype_vec] = ...
            %                MixtureOfGaussiansSimulateData([(1-f_vec)^2, 2*f_vec*(1-f_vec), f_vec^2], ...
            %                [-beta 0 beta], [1 1 1], n_samples_vec); % simulate genotype (allelic) and phenotype (QTL)
            %            cur_samples_genotype_vec = cur_samples_genotype_vec-2;
            %            beta_estimate = (cur_samples_genotype_vec*cur_samples_genotype_vec')^(-1)* ...
            %                cur_samples_genotype_vec*cur_samples_QTL_vec';
            %            stat_vec(j) = beta_estimate * sum(cur_samples_genotype_vec .* (2.*cur_samples_QTL_vec - ...
            %                beta_estimate .* cur_samples_genotype_vec)); % This is the chi-square test statistic
            
            %            beta_estimate = beta(j);
            %            f_vec_estimate = cur_contigency_table_vec(2);
            %            stat_vec(j) = beta_estimate * (2 * genotype_phenotype_prod  - ...
            %                beta_esimate * genotype_sqr);
            corr_vec(j) = beta(j) .* sqrt(2.*f_vec(j) .* (1-f_vec(j))) / sqrt(V(j));
            stat_vec(j) = corr_vec(j).*sqrt((n_samples_vec(j)-2)./(1-corr_vec(j).^2)); % t-statistic
            %            [corr_vec(j) p_vals_vec(j)] = corr(vec2column(cur_samples_QTL_vec), ...
            %                vec2column(cur_samples_genotype_vec)); % compute p-val of Pearson correlation
            p_vals_vec(j) = pvalPearson('b', corr_vec(j), n_samples_vec(j)); %  2*tcdf(-stat_vec(j),n_samples_vec(j)-2)
            
        end
        % p_vals_vec = 1-chi2cdf(stat_vec, 1); % chi-square one degree of freedom (not used. We use t-test!)
        
end % switch test stat used

s = warning('on', 'all');

debug_power=0;
if(debug_power)
    if(~isempty(optional_inds))
    observed_table = mean(observed_table) / n_samples_vec(1);
    expected_table = mean(expected_table) / n_samples_vec(1);
    
    NCP_empirical = n_samples_vec(1) * ...
        sum(sum((expected_table-observed_table).^2 ./ ...
        observed_table));
    NCP_empirical2 = n_samples_vec(1) * ...
        sum(sum(expected_table - 2*observed_table + observed_table.^2./expected_table))
    NCP_empirical3 = n_samples_vec(1) * ...
        sum(sum((expected_table-observed_table).^2 ./ ...
        expected_table));

    x_vec = NCP_empirical + ( -50:0.01:50);
    first_ind = find(x_vec>=0, 1); x_vec = x_vec(first_ind:end);
    figure; hold on;
    plot(x_vec, ncx2pdf(x_vec, 1, NCP_empirical));
    [hhh bins_loc] = hist_density(stat_vec, 50);
    plot(bins_loc, hhh, 'r');    
    plot(x_vec, ncx2pdf(x_vec, 1, NCP_empirical2), 'g');

    x_alpha = chi2inv(1-alpha_vec(1), 1);
    line([x_alpha x_alpha], [0 max(hhh)], 'color', 'k', 'linewidth', 3);
    NCP_power = 1-ncx2cdf(x_alpha, 1, NCP_empirical2)
    Empirical_power = mean(stat_vec > x_alpha);
    title(['Empirical NCP=' num2str(NCP_empirical2, 3) ...
        ', NCP-normalized=' num2str(NCP_empirical2/n_samples_vec(1), 3) ',  NCP-power=' ...
        num2str(NCP_power*100,3) '%, Empirical-power=' num2str(Empirical_power*100,3) '%']);
    
    
%    p_vals_vec(j) = 1-chi2cdf(stat_vec(j),1);
    end

end
debug_power=1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  plot_stat_figures_internal(stat_vec, z_expected, z_expected_MLT, contigency_table_vec, ...
    n_samples_vec, mu, h_x, test_stat, test_params)

figure; hold on;
hist_density(stat_vec, 100); % plot data
x_vec = linspace(min(min(stat_vec)-30, -5), max(max(stat_vec)+30, 5), 500);
%x_vec = -5:0.01:5; % :1:n_samples_vec(1);
switch test_stat
    case 'gof'
        analytic_str = 'Guass. Approx.';
        plot(x_vec, normpdf(x_vec), 'r'); % , stat_mu, stat_sigma), 'r');
    case 'gof2d'
        analytic_str = '\chi^2() Approx.';
        dof = test_params(1); non_centrality_parameter = test_params(2); % plot also the non-null distribution
        plot(x_vec, chi2pdf(x_vec, dof), 'r');
        plot(x_vec, ncx2pdf(x_vec, dof, non_centrality_parameter), 'g');
        
    case 'LRT'
        analytic_str = 'LRT Guass. Approx.';
        LRT_mu = test_params(1); LRT_var = test_params(2);
        plot(x_vec, normpdf((x_vec-n_samples_vec(1)*LRT_mu)./ ...
            sqrt(n_samples_vec(1)*LRT_var)), 'r'); % , stat_mu, stat_sigma), 'r');
        
    otherwise % other tests not implemented yet
        analytic_str = '';
        
end
xlabel([test_stat '-test-statistic']); ylabel('freq.');
title('Fit of LT model');
legend('Observed', analytic_str);   %            ['\chi^2_{' num2str(n_samples_vec(1)-1) '}']);

figure; hold on; hist_density(z_expected(:), 500);
title(['Expected Z for different genotypes under LT. Mean:' ...
    num2str(100*mean(z_expected(:)),3) ...
    '%, std:' num2str(100*std(z_expected(:)),3) '%']);
xlabel('Pr(Z=1)'); ylabel('Freq.');

z_observed = contigency_table_vec(:,:,3); z_observed = z_observed(:);
z_expected = z_expected(:); z_expected_MLT = z_expected_MLT(:);

%                figure; plot(z_expected, z_observed, '.');

num_bins = 100; %bins = linspace(0,1,num_bins+1)
bins = quantile(z_expected(:), linspace(0,1,num_bins+1));
empirical_z_observed = zeros(num_bins,1);
empirical_z_expected = zeros(num_bins,1);
empirical_z_expected_MLT = zeros(num_bins,1);
for j=1:num_bins
    good_inds = find( (z_expected >= bins(j)) & (z_expected <= bins(j+1)) );
    empirical_z_observed(j) = mean(z_observed(good_inds));
    empirical_z_expected(j) = mean(z_expected(good_inds));
    empirical_z_expected_MLT(j) = mean(z_expected_MLT(good_inds));
end
figure; hold on; % plot expected vs. observed
plot(empirical_z_expected, empirical_z_observed, '.');
plot(empirical_z_expected_MLT, empirical_z_observed, 'g.');
plot(0:1, 0:1, 'r');
legend({'LT expected', 'MLT expected', 'y=x'}, 4);
xlabel('expected z (MLT/LT model)'); ylabel('observed z  - (LT/MLT)');
title(['Data generated under LT model with \mu=' num2str(mu*100,3) ...
    '%. h_x=' num2str(100*mean(h_x),3) '% (2 loci). (' ...
    num2str(round(length(z_expected(:))/num_bins)) ' points-per-bin). n=' ...
    num2str(n_samples_vec(1))]);

X1 = mat2vec(contigency_table_vec(:,:,1));
X2 = mat2vec(contigency_table_vec(:,:,2));
num_bins2d = 10; %bins = linspace(0,1,num_bins+1) % divide to 2d bins
bins2dx1 = quantile(X1, linspace(0,1,num_bins2d+1));
bins2dx2 = quantile(X2, linspace(0,1,num_bins2d+1));
empirical_z_observed = zeros(num_bins2d);
empirical_z_expected = zeros(num_bins2d);
empirical_z_expected_MLT = zeros(num_bins2d);

for j=1:num_bins2d % loop on 2d bins
    for k=1:num_bins2d
        good_inds = find( (X1 >= bins2dx1(j)) & (X1 <= bins2dx1(j+1)) );
        good_inds = intersect(good_inds, ...
            find( (X2 >= bins2dx2(k)) & (X2 <= bins2dx2(k+1)) ));
        empirical_z_observed(j,k) = mean(z_observed(good_inds));
        empirical_z_expected(j,k) = mean(z_expected(good_inds));
        empirical_z_expected_MLT(j,k) = mean(z_expected_MLT(good_inds));
        mean_X1(j,k) = mean(X1(good_inds));
        mean_X2(j,k) = mean(X2(good_inds));
    end
end
figure; hold on; % plot expected vs. observed
plot(empirical_z_expected(:), empirical_z_observed(:), 'r.');
plot(empirical_z_expected_MLT(:), empirical_z_observed(:), 'm.');
plot(0:1, 0:1, 'r');
legend({'LT expected', 'MLT expected', 'y=x'}, 4);
xlabel('expected z (MLT/LT model)'); ylabel('observed z  - (LT/MLT)');
title(['Data generated under LT model with \mu=' num2str(mu*100,3) ...
    '%. h_x=' num2str(100*mean(h_x),3) '% (2 loci). two-dim bins. n=' ...
    num2str(n_samples_vec(1))]);
figure; hold on;
plot(mean_X1(:) + mean_X2(:), empirical_z_expected(:), 'r.');
plot(mean_X1(:) + mean_X2(:), empirical_z_expected_MLT(:), 'm.');
legend({'LT expected', 'MLT expected'}, 4);
xlabel('X1+X2'); ylabel('Z expected');
title('X1+X2 vs. z expected under two models');


figure; hold on; % plot expected z under two models
%plot(z_expected, z_expected_MLT, 'm.');
hist2d_draw([z_expected z_expected_MLT], 0:0.02:1, 0:0.02:1, ...
    'expected z (LT)', 'expected z (MLT)', ...
    'Expeted z for two alternate models');
%                 xlabel('expected z (LT)'); ylabel('expected z (MLT)');
%                 title('Expeted z for two alternate models');

figure; hold on; % Plot distribution of difference in z expected
hist_density(z_expected - z_expected_MLT, 500);
title(['Expectation diff. under two models. Mean:' ...
    num2str(100*mean(z_expected - z_expected_MLT),3) ...
    '%, std:' num2str(100*std(z_expected - z_expected_MLT),3) '%']);
xlabel('Pr(z_{LT}=1)-Pr(z_{MLT}=1)'); ylabel('freq.');
