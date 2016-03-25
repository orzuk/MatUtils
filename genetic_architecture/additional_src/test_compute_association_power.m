% Run tests for power calculations

AssignGeneralConstants;
test_epistasis = 0; % test pairwise epistasis using a simple 2x2 test (allelic)
test_power_analytic = 1; test_lambda=1; % New: test that analytic and sampling-based power calculations match
test_logistic_regression = 0; % test (deviation from?) logistic regression model
test_RVAS_power_for_paper = 1; % test rvas power calculator (program should accompany PNAS RVAS paper)


%n_samples_vec = 50:50:10000; % number of cases and controls
n_samples_vec = ceil(logspace(1, 6, 100));
alpha = 0.05/20000; % 0.0005; % 5*10^(-8); % 0.000005; % confidence level required
iters = 1000;
theta = 0.5; % proportion of cases in sample


%% NEW! Test RVAS power 
if(test_RVAS_power_for_paper) 
    power_parameters_file = 'C:\research\matlab\libs\genetic_architecture\Data\powerParameters.txt';
    power_output  = RVAS_power_calculator(power_parameters_file);
    
end
    

p_mat = cell(20,1);
if(test_lambda)
%    lambda_vec = [2 5]-1;
    lambda_vec = [1.5 2 3 5 10 20]-1; % set different relative-risk values. Only strong effects !!!
    f_rare = 0.001; % take rare alleles
    
    
    for i=1:length(lambda_vec)
        p_mat{i} = genetic_relative_risk_to_p_z_x_marginal(f_rare, lambda_vec(i)+1, 0.05);
    end
    
else
    p_mat{1} = [0.81 0.09; ...
        0.09 0.01]; % no effect. Power should be equal to alpha
    p_mat{2} = [0.81 0.085; ...
        0.085 0.02]; % moderate effect. Power should be bigger than alpha
    p_mat{3} = [mat_into_vec(p_mat{1}) mat_into_vec(p_mat{1})]; % large effect, but no epistatic effect. Power should be big for single-locus but ~alpha for epistatic test
    
    % Switch 0 and 1 - shouldn't change results!
    p_mat{2} = p_mat{2}(:,2:-1:1); % flipped matrix with effect size zero !!!
    
    
    % New: load a real SNP's statistics
    %%%%%%%load('c:\public_html\data\common_disease_broad_catalog\disease_data.mat');
    % load(fullfile(html_outdir, 'common_disease_broad_catalog\disease_data.mat'));
    %%%%%%%p_mat{2} = vec2mat(genetic_relative_risk_to_p_z_x_marginal(data.MAF(577), ...
    %%%%%%%    data.OR(577), data.Prevalence(577)), 2);
    
    %p_mat{2} = [0.05 0.85; 0.005 0.085] ./ 0.99; % similar matrix with effect size zero !!!
    
    
    
    % p_mat{2} = vec2mat(genetic_relative_risk_to_p_z_x_marginal(0.9, 1, 0.1), 2) % put no effect, high f_vec (e.g. 0.9) to get a bad test! (high power : false positives)
    %p_mat{2} = vec2mat(genetic_relative_risk_to_p_z_x_marginal(0.1, 1.4, 0.1), 2) % put high f_vec (e.g. 0.9) to get a bad test! (high power : false positives)
    % p_mat{2} = [0.54 0.06; 0.36 0.04]; % no effect. High MAF but still lower than 0.5
    
    model_type = 'logistic'; % 'additive' %
    test_stat = 'chi-square';
    switch model_type
        case 'additive'
            p_z_given_x = [0.1 0.12; 0.13 0.15]; % additive model, small log ratio
        case 'multiplicative'
            p_z_given_x = [0.1 0.2; 0.3 0.6]; % multiplicative model
            %        p_z_given_x = [0.1 0.12; 0.13 0.16]; % additive model, small log ratio
            
        case 'logistic'
            beta = 0.5; alpha_1 = 0.523; alpha_2 = -1.213;
            p_z_given_x = exp( beta + alpha_1 .* [0 0; 1 1] + alpha_2 .* [0 1; 0 1]);
            p_z_given_x = p_z_given_x ./ (1+p_z_given_x)
    end
    p_mat{3} = p_mat{3} .* [1-mat_into_vec(p_z_given_x) mat_into_vec(p_z_given_x)];
    
end % if use effect size lambda



if(test_epistasis)
    for i=4:20
        p_mat{i} = [mat_into_vec(p_mat{1}) mat_into_vec(p_mat{1})]; % large effect, also epistatic effect. Power should be big for single-locus but grow slower for epistatic test
        p_z_given_x = rand(2); % [0.1 0.2; 0.3 0.6]; % strong multiplicative model
        p_z_given_x = p_z_given_x ./ sum(p_z_given_x(:));
        p_mat{i} = p_mat{i} .* [1-mat_into_vec(p_z_given_x) mat_into_vec(p_z_given_x)];
        p_vec{i} = vec2row(mat_into_vec(vec2row(p_mat{i}))); % vector of 1x4 or 1x8 probabilities
        if(length(p_vec{i}) == 8) % check if indeed the model is epistatic or not
            [epistasis_flag_original epistasis_deviation(i)] = ...
                test_model_for_epistasis(p_vec{i}, test_stat, epsilon);
            [epistasis_flag_case_control epistasis_deviation_case_control(i)] = ...
                test_model_for_epistasis(pop_prob_to_case_control_prob(p_vec{i}), ...
                test_stat, epsilon);
        end
        
    end
    
    
    figure; plot(epistasis_deviation, epistasis_deviation_case_control, '.');
    hold on; plot([min(epistasis_deviation):0.001:max(epistasis_deviation)], ...
        [min(epistasis_deviation):0.001:max(epistasis_deviation)], 'r');
    xlabel('orig. deviation'); ylabel('case-control deviation');
    title('epistasis effect before/after case-control transformation');
    signal_type = {'no-effect', 'small-marginal-large-epistasis', ...
        'only-marginal-no-epistatis', 'large-marginal-small-epistasis'};
    sampling_type = {'all_population', 'case-control'}
    
    test_stat =  'chi-square'; % 'hypergeometric'; % 'chi-square'; % 'hypergeometric'; % 'chi-square'; % 'hypergeometric'; % 'chi-square', 'z-score'
    test_type = 'single-locus'; % 'epistasis';
    pow = cell(3,1);
    test_different_powers = 0;
    if(test_different_powers)
        for i=3:length(p_mat)
            if(length(p_mat{i}) == 2)
                test_type = 'single-locus';
                test_stat = 'chi-square';
            else
                test_type = 'epistasis';
                %       test_stat = 'logistic';
                test_stat = 'chi-square';
            end
            figure; hold on;
            color_ctr=1;
            for k=1:2 % case-control or random population
                for j=1:2 % do single locus and epistasis
                    if(j == 1)
                        test_p_mat = p_mat{i}(1:2,:) + p_mat{i}(3:4,:); % collapse model to single locus
                        test_type = 'single-locus';
                        test_stat = 'chi-square';
                    else
                        test_type = 'epistasis';
                        test_stat = 'logistic'; % 'chi-square'; % 'logistic';
                        test_p_mat = p_mat{i};
                    end
                    
                    pow{i,j} = compute_association_power(test_p_mat, n_samples_vec, [], alpha, ...
                        iters, test_type, test_stat, sampling_type{k});
                    plot(n_samples_vec, pow{i,j}, color_vec(color_ctr), 'linewidth', 2);
                    plot(n_samples_vec, pow{i,j}, [color_vec(color_ctr) 'x']);
                    color_ctr=color_ctr+1;
                end % loop on j
            end
            if(i <= 4)
                title(['power to detect effect for various sample sizes at \alpha=' ...
                    num2str(alpha) ' for ' signal_type{i} ' effect']);
                xlabel('N'); ylabel('power');
                if(color_ctr==5)
                    legend({'main-effect-population', '', 'epistasis-population', '', ...
                        'main-effect-case-control', '', 'epistasis-case-control', ''}, ...
                        'location', 'northwest');
                else
                    legend({'main-effect', '', 'epistasis', ''}, 'location', 'northwest');
                end
            end
        end
    end % test all different power calculations
    power_is = pow
    
end % test epistasis


if(test_power_analytic)
   
    for i=1:length(lambda_vec) %  2:2 % non-epistatic model 3:length(p_mat)
        num_cases = round(n_samples_vec*theta); num_controls = n_samples_vec - num_cases;
        test_p_mat = p_mat{i}; % (1:2,:) + p_mat{i}(3:4,:); % collapse model to single locus
        if(~test_lambda)
            test_p_mat = test_p_mat(:,2:-1:1);
        end
        [power_sampling{i} pvals_vec{i} chi_stat_vec{i}] = ...
            compute_association_power(test_p_mat, num_cases, num_controls, alpha, ...
            iters, 'single-locus', 'chi-square', 'case-controls');
        [power_analytic{i}, ~, ~, non_centrality_parameter{i}] = ...
            compute_association_power(test_p_mat, num_cases, num_controls, alpha, ...
            iters, 'single-locus', 'chi-square-analytic', 'case-controls');
        
        if(9999 < 0) % additional tests (exclude for now)
            [power_Armitage_sampling{i} pvals_Armitage_vec{i} chi_stat_Armitage_vec{i}] = ...
                compute_association_power(test_p_mat, num_cases, num_controls, alpha, ...
                iters, 'armitage', 'chi-square', 'case-controls');
            [power_Armitage_analytic{i}, ~, ~, non_centrality_Armitage_parameter{i}] = ...
                compute_association_power(test_p_mat, num_cases, num_controls, alpha, ...
                iters, 'armitage', 'chi-square-analytic', 'case-controls');
            [power_Armitage_analytic2{i}, ~, ~, non_centrality_Armitage_eliana_parameter{i}] = ...
                compute_association_power(test_p_mat, num_cases, num_controls, alpha, ...
                iters, 'armitage-eliana', 'chi-square-analytic', 'case-controls');
        end
        
        %       figure; % NEW! Plot power valuse for case-only study
        
        
        %        break;
        
        if(9999 < 0)
            
            figure; hold on; % Plot power values
            plot(n_samples_vec, power_sampling{i}, 'xb', 'linewidth', 2);
            plot(n_samples_vec, power_Armitage_sampling{i}, 'xg', 'linewidth', 2);
            plot(n_samples_vec, power_analytic{i}, 'rx');
            plot(n_samples_vec, power_Armitage_analytic{i}, 'xc', 'linewidth', 2);
            legend('sampling', 'Armitage-sampling', 'analytic', 'Armitage-analytic');
            title('Power calculated for different sample sizes');
            xlabel('num samples'); ylabel('power');
            
            figure; hold on;
            plot(power_analytic{i}, power_sampling{i}, 'xb', 'linewidth', 2);
            plot(power_Armitage_analytic{i}, power_Armitage_sampling{i}, 'xg', 'linewidth', 2);
            plot(0:0.01:1, 0:0.01:1, 'r');
            title('Power calculated for different sample sizes');
            legend({'binary', 'armitage', 'Y=X'}, 4);
            xlabel('power (analytic)'); ylabel('power (sampling)');
            
            figure; hold on; % Plot non-central chi-square distribution
            hist_density(chi_stat_vec{i}, 100, 'r', 1, 1, 1); % Plot empirical distribution
            hist_density(chi_stat_Armitage_vec{i}, 100, 'g', 1, 1, 1); % Plot empirical distribution armitage
            x_vec = 0:0.01:max(chi_stat_vec{i}) * 3.1;
            plot(x_vec,  ncx2pdf(x_vec, 1, non_centrality_parameter{i})); % plot non centralized chi-square distribution
            plot(x_vec,  ncx2pdf(x_vec, 1, non_centrality_Armitage_parameter{i}), 'm'); % plot non centralized chi-square distribution
            plot(x_vec,  ncx2pdf(x_vec, 1, non_centrality_Armitage_eliana_parameter{i}), 'c'); % plot non centralized chi-square distribution
            title('Empirical vs. Theoretical Statistic Distribution');
            xlabel('Val'); ylabel('Freq.'); legend('empirical', 'empirical-armitage', ...
                'theoretical', 'theoretical-armitage', 'theoretical-armitage-eliana');
        end % if 999<0
    end % loop on models
    
    figure; % NEW! Plot power valuse for case-only study
power_analytic_mat = cell2vec(power_analytic);
power_sampling_mat = cell2vec(power_sampling);
%%%    for i=1:length(lambda_vec) %  2:2 % non-epistatic model 3:length(p_mat)
h1 =         semilogx(n_samples_vec, power_analytic_mat, 'linewidth', 2); hold on;
%%%    end
%%%    for i=1:length(lambda_vec) %  2:2 % non-epistatic model 3:length(p_mat)
        semilogx(n_samples_vec, power_sampling_mat, '--', 'linewidth', 2); % hold on;
        
%%%    end
    %                            pos_l = get(h_leg, 'position'); set(h_leg, 'position', [pos_l(1)+0.2835 pos_l(2)-0.1435 pos_l(3) pos_l(4)]);
        c1 = get(h1, 'color');

    for i=1:length(lambda_vec) % TEXT 
        tmp_ind =  find(power_analytic{i}>0.5, 1); 
        if(isempty(tmp_ind))
            tmp_ind=1;
        end
        N50(i) = n_samples_vec(tmp_ind);
        text(1+N50(i)*0.5, 0.5, num2str(1+lambda_vec(i)), 'color', c1{i}, 'fontsize', 14); 
    end
    text(1+N50(end)*0.25, 0.5, '1+\lambda=', 'color', 'k', 'fontsize', 14); 
    
        
%    legend_vec = num2str(lambda_vec');
%    h_leg = legend(legend_vec, 2); % {'Sampling', 'Analytic'},
%    set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
%    title('Power calculated for different sample sizes');
    xlabel('Num cases'); ylabel('Power');
    xlim([10^2 10^6]);
    add_faint_grid(0.5);
    power_figs_dir = '../../common_disease_model/figs/EyreWalker/new_eric/power';
    my_saveas(gcf, fullfile(power_figs_dir, 'power', 'NCP_approx'), {'epsc', 'pdf', 'jpg'}); % save figure
    
end % test power analytic vs. sampling


if(test_logistic_regression)
    N=100;
    beta = 0.5; alpha_1 = 0; alpha_2 = -2.213;
    X = rand(N,2) > 0.8; % test logistic regression
    Y = ilogit(beta + alpha_1 .* X(:,1) + alpha_2 .* X(:,2));
    Y = rand(N,1) < Y;
    
    [b dev stats] =  glmfit(X, [Y ones(N,1)], 'binomial', 'link', 'logit');
end




test_ncp = 1; % test computation lf non-centrality parameter
if(test_ncp)
    grr = 1.2; n_cases = 2000; n_controls = 4000;
    ncp_standard = 2*genetic_relative_risk_to_non_centrality_parameter(0.1, grr, n_cases, n_controls, 0.1, 'standard_chi_square')
    ncp_apples = genetic_relative_risk_to_non_centrality_parameter(0.1, grr, n_cases, n_controls, 0.1, 'apples')
    ncp_eliana = genetic_relative_risk_to_non_centrality_parameter(0.1, grr, n_cases, n_controls, 0.1, 'eliana')
    ncp_ott = genetic_relative_risk_to_non_centrality_parameter(0.1, grr, n_cases, n_controls, 0.1, 'ott')
end


test_qtl = 0; % test analytical and empirical power for QTL
if(test_qtl)
    for    beta = [0 0.1]
        iters = 1000; % heavy part ..
        RAF = 0.2;
        V = 1;
        alpha = 0.05
        test_p_mat_QTL = [RAF beta 0 V];
        n_samples_vec = 100:100:500;
        [power_QTL_sampling pvals_vec_QTL chi_stat_vec_QTL] = ...
            compute_association_power(test_p_mat_QTL, n_samples_vec, [], alpha, ...
            iters, 'single-locus', 'chi-square-QTL', 'population');
        [power_QTL_analytic, ~, ~, non_centrality_parameter_QTL] = ...
            compute_association_power(test_p_mat_QTL, n_samples_vec, [], alpha, ...
            iters, 'single-locus', 'chi-square-QTL-analytic', 'population');
        figure; hold on; plot(power_QTL_analytic, power_QTL_sampling, '.');
        plot(0:0.01:1, 0:0.01:1, 'r');
        title(['QTL linear test power (analytic vs. sampling \beta = ' num2str(beta)]);
        xlabel('Analytic power'); ylabel('Sampling power');
        
        figure; hold on; % Plot non-central chi-square distribution
        hist_density(chi_stat_vec_QTL, 100, 'r', 1, 1, 1); % Plot empirical distribution
        x_vec = (min(chi_stat_vec_QTL):0.01:max(chi_stat_vec_QTL)) * 1.1;
        plot(x_vec,  nctpdf(x_vec, 1, non_centrality_parameter_QTL(end))); % plot non centralized chi-square distribution
        title(['Empirical vs. Theoretical QTL Statistic Distribution \beta = ' num2str(beta)]);
        xlabel('Val'); ylabel('Freq.'); legend('empirical', 'theoretical');
    end % loop on beta
end
