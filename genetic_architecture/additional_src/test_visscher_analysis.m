% Test the sibling regression method from Visscher et al. (2006) for the LP model
%function test_visscher_analysis()
figs_dir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/';
IBD_file = '../../common_disease_model/data/visscher/IBD_data.mat';
IBD_dir = '../../common_disease_model/data/visscher';

run_full=0; % compute r_R for full range of IBD sharing
run_power=0; % compute power to detect slope different from epidemiological estimate
run_one=0; % run only around half (siblings)
run_binary=0; % run for disease
run_shared_env=0; % here run models such that at 1/2 they all have the same shared environment
run_schizophrenia_data=1; % fit schizophrenia example 
debug_MZ=0;
h_x = 0.8; % (narrow-sense) heritability
r_DZ = h_x/2;
mu=0.01; % prevalence for binary trait
iters = 1000000;
n_samples_vec = round(logspace(3,6,20)); %   [20000:20000:100000];
%n_samples_vec = 1000:1000:10000; % [10000:10000:100000];
IBD_mean = 0.5;
IBD_std = 0.036; % value reported by Visscher for sibs % 0.1;
lambda_mz = heritability_to_familial_risk(h_x, 'liability', mu, 1);

ttt=cputime;
if(run_one) % run around half (siblings) 
    %figure; subplot(2,1,1); hold on;
    [BETA h_loci] = simulate_visscher_analysis_LP(iters, 1, h_x, [], IBD_mean, IBD_std); % additive model
    [BETA_again, IBD_sharing_vec, qtl_R, ~,power_vec] = ... % compute power !!!
        simulate_visscher_analysis_LP(n_samples_vec, 1, h_x, [], IBD_mean, IBD_std,[],'sampling') % additive model
    figure; plot(n_samples_vec, power_vec);
    title('Power to distinguish model from additivity using IBD slope for siblings');
    xlabel('sample size'); ylabel('power');
    
    
    %subplot(2,1,2); hold on;
    
    %figure;
    [BETA3 h_loci3] = simulate_visscher_analysis_LP(iters, ...
        3, h_x, [], IBD_mean, IBD_std); % LP(3) model
    [BETA3_again] = simulate_visscher_analysis_LP(iters, ...
        3, h_x, [], IBD_mean, IBD_std,[],'sampling') % LP(3) model
    [BETA2_again] = simulate_visscher_analysis_LP(iters, ...
        2, h_x, [], IBD_mean, IBD_std,[],'sampling') % LP(3) model
    
    %my_saveas(gcf, ...
    %    fullfile(figs_dir, 'supp_fig7visscher_analysis'), 'epsc');
    %close all
    
    
end % run one (around siblings) 

N_vec = [1:5 10]; num_N=length(N_vec);
power_vec = zeros(length(N_vec),length(n_samples_vec)); power_vec_numeric=power_vec;
h_x_one_liab = zeros(length(N_vec),1); h_x_one_liab_disease = h_x_one_liab;
for i=1:length(N_vec) % loop on N 
    N=N_vec(i)
    run_N=i
    if((run_full ||  run_shared_env) || (run_power || run_binary)) % fitting parameters (heavy computation)
        if((~exist('h_x_one_liab', 'var') || (length(h_x_one_liab) < i)) ...
                || (h_x_one_liab(i)==0))
            h_x_one_liab(i) = familial_risk_to_LP_parameters(h_x, N,  1, 'quantitative');
        end
        if((~exist('h_x_one_liab_shared', 'var') || (length(h_x_one_liab_shared) < i)) ...
                || (h_x_one_liab_shared(i)==0))
            [h_x_one_liab_shared(i) c_R_shared(i)] = ...
                familial_risk_to_LP_parameters([h_x r_DZ], N,  [1 0.5], 'quantitative');
        end
    end
    if(run_full)
        [BETA_full(i), IBD_sharing_vec(i,:), qtl_R(i,:), ~, ~, ...
            h_loci_full(i), h_pop(i)] = ...
            simulate_visscher_analysis_LP(10, N, h_x_one_liab(i), [], ...
            [], [], 'full-range'); % make the MZ correlation the same
    end
    if(run_shared_env) % run model with shared environment set such that r_MZ and r_DZ remain fixed
        [BETA_shared(i), IBD_sharing_vec_shared(i,:), qtl_R_shared(i,:), ~, ~, ...
            h_loci_shared(i), h_pop_shared(i)] = ...
            simulate_visscher_analysis_LP(10, N, h_x_one_liab_shared(i), c_R_shared(i), ...
            [], [], 'sibling', 'numeric'); % make the MZ correlation the same
    end
    if(run_power)
        %         [~, ~, ~, power_vec(i,:)] = ...
        %             simulate_visscher_analysis_LP(n_samples_vec, N, h_x_one_liab(i), [], ...
        %             [], [], 'sibling', 'sampling'); % Just compute power
        [~, ~, ~, power_vec_numeric(i,:)] = ...
            simulate_visscher_analysis_LP(n_samples_vec, N, h_x_one_liab(i), [], ...
            IBD_mean, IBD_std, 'sibling', 'numeric'); % Just compute power
    end
    if(run_binary)
        if((~exist('h_x_one_liab_disease', 'var') || (length(h_x_one_liab_disease) < i)) ...
                || (h_x_one_liab_disease(i)==0))
            h_x_one_liab_disease(i) = familial_risk_to_LP_parameters(lambda_mz, N,  1, 'binary', mu);
        end
        [BETA_binary(i), IBD_sharing_vec(i,:), qtl_R_binary(i,:), lambda_R_binary(i,:), ~,...
            h_loci_binary(i), h_pop_binary(i)] = ...
            simulate_visscher_analysis_LP(10, N, h_x_one_liab_disease(i), [], ...
            [], [], 'full-range', 'numeric', 'binary', mu); % make the MZ correlation the same
    end
end % loop on N 
%save(IBD_file, '
%load(IBD_file);

%mid_ind = ceil(0.5*size(qtl_R, 2))
if(exist('BETA_full', 'var')) % run 'full' model -  IBD sharing also around zero
    figure; % plot again
    ylim([0 h_x*1.001]);
    title(['IBD sharing vs. corr(Z) for LP(k,r_{MZ}=' num2str(100*h_x,3) '%)']);
    xlabel('IBD'); % ylabel('corr(Z)');
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig8IBD_vs_phenotype_correlation_empty'), {'epsc', 'pdf'});
    plot(IBD_sharing_vec', qtl_R', 'linewidth', 2);
    legend([repmat('k=', length(N_vec), 1) num2str((N_vec)')],4);
    ylim([0 h_x*1.001]);
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig8IBD_vs_phenotype_correlation'), {'epsc', 'pdf'});
    %close all
    
    
    figure; hold on;
    plot(100*BETA_full, 100*min(1,h_pop), '.'); % Print slope at one half
    xlabel('h_{IBD} (%)'); ylabel('h_{pop} (ACE)  (%)'); % qtl_R(:,mid_ind)
    for i=1:length(N_vec)
        text(100*BETA_full(i)+0.1, 100*min(1,h_pop(i))+0.1, ['N=' num2str(N_vec(i))]);
    end
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig9h_IBD_vs_h_pop'), 'epsc');
    
    figure; hold on;
    plot(N_vec, 100*h_loci_full, '*b'); % Print h_loci
    plot(N_vec, 100*BETA_full, '*r'); % Print slope at one half (h_IBD)
    plot(N_vec, 100*min(1,h_pop), '*g'); % Print h_pop
    legend('h_{all}', 'h_{IBD}', 'h_{pop}');
    plot(N_vec, 100*h_loci_full, 'b'); % Print h_loci
    plot(N_vec, 100*BETA_full, 'r'); % Print slope at one half (h_IBD)
    plot(N_vec, 100*min(1,h_pop), 'g'); % Print h_pop
    xlabel('N'); ylabel('h estimator  (%)'); % qtl_R(:,mid_ind)
    for i=1:length(N_vec)
        text(100*BETA_full(i)+0.1, 100*min(1,h_pop(i))+0.1, ['N=' num2str(N_vec(i))]);
    end
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig10h_loci_vs_h_IBD_vs_h_pop'), 'epsc');
    
end % plot full model

if(run_shared_env) % Plot model with shared environment (IBD variation just around 1/2)
    figure; hold on;
    plot(N_vec, 100*h_loci_shared, '*b'); % Print h_loci
    plot(N_vec, 100*BETA_shared, '*r'); % Print slope at one half (h_IBD)
    plot(N_vec, 100*min(1,h_pop_shared), '*g'); % Print h_pop
    legend('h_{all}', 'h_{dr_s}', 'h_{pop}', 3);
    plot(N_vec, 100*h_loci_shared, 'b'); % Print h_loci
    plot(N_vec, 100*BETA_shared, 'r'); % Print slope at one half (h_IBD)
    plot(N_vec, 100*min(1,h_pop_shared), 'g'); % Print h_pop
    xlabel('N'); ylabel('h estimator  (%)'); % qtl_R(:,mid_ind)
    for i=1:length(N_vec)
        text(100*BETA_shared(i)+0.1, 100*min(1,h_pop_shared(i))+0.1, ['N=' num2str(N_vec(i))]);
    end
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig12h_loci_vs_h_IBD_vs_h_pop_with_shared_env'), 'epsc');
    
    R = [vec2column(N_vec) 100*[repmat(h_x, num_N, 1) repmat(r_DZ, num_N, 1) ...
        vec2column(h_loci_shared)  vec2column(h_pop_shared) vec2column(BETA_shared) ...
        vec2column(h_x_one_liab_shared) vec2column(c_R_shared)]];
    R(R<0.01)=0; % get rid of numerical erros
    R = [{'N', 'r_{MZ} (%)', 'r_{DZ} (%)', ...
        'h_{all} (%)', 'h_{pop} (%)',  'h_{dr_s}', ...
        'h_{pathway} (%)', 'c_{R} (%)'}' ...
        num2str_cell(num2cell(R), 3)']';
    shared_file = fullfile(IBD_dir, 'quant_trait_shared_env.txt');
    savecellfile(R, shared_file,[],1); % save value for shared environment
end



if(run_power) % Compute power % exist('power_vec', 'var'))
    figure; % hold on;
    %    plot(n_samples_vec, power_vec); hold on; % semilogx
    semilogx(n_samples_vec, power_vec_numeric); % , '--'); % semilogx
    legend([repmat('N=', length(N_vec), 1) num2str((N_vec)')],2);
    xlabel('Sample size'); ylabel('Power');
    %    legend('sampling', 'numeric');
    title('Power to detect difference between h_{IBD} (slope) and h_{pop} for LP(N,h_x) model');
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig11power_h_IBD_vs_h_pop'), 'epsc');
end

if(exist('BETA_binary', 'var')) % plot for binary variables
    figure; % plot lambda_R vs. IBD sharing
    plot(IBD_sharing_vec', lambda_R_binary');
    legend([repmat('N=', length(N_vec), 1) num2str((N_vec)')],4);
    %    ylim([0 h_x*1.001]);
    title(['IBD sharing vs. \lambda_R for LP(N,r_{MZ}=' num2str(100*h_x,3) '%, \mu=' num2str(mu) ')']);
    xlabel('IBD'); ylabel('\lambda_R');
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig11IBD_vs_lambda_R'), 'epsc');
    figure; % plot inferred qtl_R vs. IBD sharing
    plot(IBD_sharing_vec', qtl_R_binary');
    legend([repmat('N=', length(N_vec), 1) num2str((N_vec)')],4);
    ylim([0 h_x*1.001]);
    title(['IBD sharing vs. inferred r_R for LP(N,r_{MZ}=' num2str(100*h_x,3) '%, \mu=' num2str(mu) ')']);
    xlabel('IBD'); ylabel('r_R');
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig11IBD_vs_qtl_R_infererd_from_lambda_R'), 'epsc');
    
    
    figure; hold on;
    plot(N_vec, 100*h_loci_binary, '*b'); % Print h_loci
    plot(N_vec, 100*BETA_binary, '*r'); % Print slope at one half (h_IBD)
    plot(N_vec, 100*min(1,h_pop_binary), '*g'); % Print h_pop
    legend('h_{all}', 'h_{dr_s}', 'h_{pop}',3);
    plot(N_vec, 100*h_loci_binary, 'b'); % Print h_loci
    plot(N_vec, 100*BETA_binary, 'r'); % Print slope at one half (h_IBD)
    plot(N_vec, 100*min(1,h_pop_binary), 'g'); % Print h_pop
    xlabel('N'); ylabel('h estimator  (%)'); % qtl_R(:,mid_ind)
    for i=1:length(N_vec)
        text(100*BETA_binary(i)+0.1, 100*min(1,h_pop_binary(i))+0.1, ['N=' num2str(N_vec(i))]);
    end
    title('Different heritability estiamtors for LP(N,h_x,\mu=0.01) disease model');
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig11h_loci_vs_h_IBD_vs_h_pop_disease'), 'epsc');
    
end



if(debug_MZ)  % Test problem with computing r_MZ (temporary)
    N=10
    qtl_familial_correlations_internal(N, N, 1, 'MAX', [], h_x_one_liab(end), 0, 'numeric', [], 'MZ')
    familial_risk_to_heritability( qtl_R(end,end), 'liability', mu, 1)
    mu_l = 1-(1-mu).^(1/N);
    
    h_mz = 0.8;
    lambda_mz = heritability_to_familial_risk(h_mz, 'liability', mu, 1);
    [lambda_R STATS] = ...             % Need to extract h_loci from here!
        compute_k_of_N_liabilities_statistics(N, 1, mu_l, h_x_one_liab(end), 0, 1);
    
    [h_x_one_liab_qtl1 c_R_qtl1] = familial_risk_to_LP_parameters(h_mz, 1,  1, 'quantitative')
    [h_x_one_liab1 c_R1] = familial_risk_to_LP_parameters(lambda_mz, 1,  1, 'binary', mu)
    
    [h_x_one_liab_qtl c_R_qtl] = familial_risk_to_LP_parameters(h_mz, N,  1, 'quantitative')
    [h_x_one_liab c_R] = familial_risk_to_LP_parameters(lambda_mz, N,  1, 'binary', mu)
end


if(run_schizophrenia_data) % Fit data from Wray&Visscher 2009 Schiz. paper
    fit_str_vec =[];
    ctr=1;
    prevalence = [0.0085 0.00407]; % get two prevalences
    lambda_mz = 52.1;
    lambda_s = [8.6; 8.55]; % 
%    lambda_s = [8.6 10.0 14.2; 9.43 10.3 8.55]; % [8.6; 8.55]; % 
    lambda_grand = [3.1 3.2 3.3 3.5 3.3; 2.52 2.71 2.95 3.04 3.8]
    lambda_cousin = [1.8 2.29];
    k_R_vec = [1 0.5 0.25 0.125]; % kinship coefficient vector
    max_family_degree=4; % how many risks to compute
    
    N_vec = 1:10; % try different LP models
    num_N = length(N_vec);
    mu = prevalence; % repmat(prevalence, num_N, 1);
    for mu_ind = 1:2 % loop on two different prevalences
        for i=1:length(N_vec)
            mu_l(mu_ind,i) = 1-(1-mu(mu_ind)).^(1/N_vec(i));
        end
        for lambda_s_ind = 1:size(lambda_s,2) % 3
            for fit_type=1:3
                switch fit_type % loop on 3 different fit types 
                    case 1 % no shared environment, just fit lambda_s
                        lambda_R = lambda_s(mu_ind,lambda_s_ind);
                        k_R = k_R_vec(2);
                        fit_str = ['mu=' num2str(mu(mu_ind)*100) '%, lambda_s=' num2str(lambda_s(mu_ind,lambda_s_ind)) ', c_R=0'];
                        c_R_vec = [];
                    case 2 % with shared environment, fit lambda_s and lambda_grand
                        lambda_R = [lambda_s(mu_ind,lambda_s_ind) lambda_grand(mu_ind,3)];
                        k_R = k_R_vec(2:3);
                        fit_str = ['mu=' num2str(mu(mu_ind)*100) '%, lambda_s=' num2str(lambda_s(mu_ind,lambda_s_ind)) ...
                            ', lambda_grand=' num2str(lambda_grand(mu_ind,3))];
                        c_R_vec = [];
                    case 3 % New: allow different shared environment for EACH relative. What are the constraints?
                        lambda_R = [lambda_s(mu_ind,lambda_s_ind) lambda_grand(mu_ind,3)];
                        k_R = k_R_vec(2:3);
                        fit_str = ['mu=' num2str(mu(mu_ind)*100) '%, lambda_s=' num2str(lambda_s(mu_ind,lambda_s_ind)) ...
                            ', lambda_grand=' num2str(lambda_grand(mu_ind,3)), ', c_R>0 for sibs-only'];
                        c_R_vec = [1 0]; % shared enivronment exists only up to sibs
                end % switch fit type
                for i=1:length(N_vec) % fit parameters
                    run_i = i
                    all_N = length(N_vec)
                    [h_pathway(ctr) c_R(ctr)] = familial_risk_to_LP_parameters( ...
                        lambda_R, N_vec(i),  k_R, 'binary', prevalence(mu_ind), c_R_vec); % fit paremters
                    h_shared_env_one_liab(ctr) = c_R(ctr) * (1-h_pathway(ctr));
                    
                    % New! Compute shared environment on the one liability scale
                    [cur_lambda_R_no_shared_env] = ...
                        compute_k_of_N_liabilities_statistics(N_vec(i), 1, mu_l(mu_ind,i), ...
                        h_pathway(ctr), 0, max_family_degree);
                    
                    compute_k_of_N = 99
                    [cur_lambda_R cur_stat] = ...
                        compute_k_of_N_liabilities_statistics(N_vec(i), 1, mu_l(mu_ind,i), ...
                        h_pathway(ctr), h_shared_env_one_liab(ctr), max_family_degree);
                    cur_liab_no_shared_env = familial_risk_to_heritability(cur_lambda_R_no_shared_env(1), 'liability', mu(mu_ind), 1);
                    cur_liab = familial_risk_to_heritability(cur_lambda_R(1), 'liability', mu(mu_ind), 1);
                    h_shared_env(ctr) = cur_liab - cur_liab_no_shared_env;
                    
                    update_parameters = 88
                    
                    %            c_R(ctr) = c_R(ctr) / (1-h_pathway(ctr)); % our parameterization: c_R is RELATIVE contribution
                    h_all(ctr) = cur_stat.h_liab_loci;
                    h_pop(ctr) = cur_stat.h_liab_twins;
                    h_phantom(ctr) = cur_stat.h_phantom_vec;
                    fit_str_vec{ctr} = fit_str;
                    all_N_vec(ctr)=N_vec(i);
                    all_mu(ctr)=prevalence(mu_ind);
                    lambda_mz_vec(ctr) = cur_lambda_R(1);
                    lambda_s_vec(ctr) = cur_lambda_R(2);
                    lambda_grand_vec(ctr) = cur_lambda_R(3);
                    lambda_cousin_vec(ctr) = cur_lambda_R(4);
                    ctr=ctr+1;
                end % loop on N
            end % loop on fit type (with/without shared environment)
        end % loop on lambda_s paper
    end % loop on prevalence paper
    fit_str_lengths = length_cell(fit_str_vec); % Pad strings
    fit_str_max_length = max(fit_str_lengths);
    for i=1:length(fit_str_vec)
        fit_str_vec{i} = [fit_str_vec{i} repmat(' ', 1, fit_str_max_length-fit_str_lengths(i))];
    end
    %    fit_str_vec = cell2mat(fit_str_vec');
    h_shared_pathway_vec = c_R .* (1-h_pathway); % absolute (rather than relative) shared environment
    h_phantom = max(h_phantom,0);
    h_phantom(h_phantom < 0.0001) = 0; % round
    S = struct('N', vec2column(all_N_vec), ...  %'fit', vec2column(fit_str_vec), ...
        'mu', vec2column(100*all_mu), 'h_pathway', vec2column(100*h_pathway), ...
        'c_R', vec2column(100*c_R), 'h_shared_pathway', vec2column(100*h_shared_env_one_liab), ...     % parameters
        'lambda_mz', lambda_mz_vec, 'lambda_s', lambda_s_vec, ...
        'lambda_grand', lambda_grand_vec, 'lambda_cousin', lambda_cousin_vec, ... % risks
        'h_all', 100*h_all, 'h_pop', 100*min(1,h_pop), 'h_phantom', 100*h_phantom, ...
        'h_shared_env', vec2column(100*h_shared_env));
    for i=1:length(fit_str_vec)
        S.fit{i} = fit_str_vec{i};
    end
    P = [length(fields(S)) 1:length(fields(S))-1]; S = orderfields(S, P);
    field_str = fields(S)
    for i=1:length(fields(S)) % transfer to strings
        eval_str = ['S.' field_str{i} ' = num2str_cell(S.' field_str{i} ', 3);']
        eval(eval_str);
    end
    % Make space every
    
    schizophrenia_outfile = fullfile(IBD_dir, 'visscher_wray_schizophrenia_new_many_N.txt');
    cell_to_mat=1;
    S_data1 = {'Risch et al.', '', '0.85', '', '', '', '52.1', '8.6', '3.3', '1.8', '', '', ''};
    S_data2 = {'Lichtenstein et al.', '', '0.41', '', '', '', '', '8.55', '2.95', '2.29', '', '', ''};
    
    S_Cell = WriteDataFile(S, schizophrenia_outfile, cell_to_mat,[], num_N); % write results to text file 
    num_lines = size(S_Cell,1);
    S_Cell = [S_Cell(1,:)' S_data1' S_Cell(2:(num_lines+1)/2,:)' S_data2' S_Cell((num_lines+3)/2:end,:)']';
    savecellfile(S_Cell, schizophrenia_outfile);
    R_header = {'$\lambda_{MZ}$', '$\lambda_{s}$', '$\lambda_{g}$', '$\lambda_{cousin}$', ...
        '$h_{all}^2$', '$h_{pop}^2$', '$h_{pathway}^2$', '$c_R$', '$\pi_{phantom}'}; % , ...
%        'Risch et al.', '\LT', '\Lpd(2)', 'Lichtenstein et al.', '\LT', '\LPd(2)'}; 
    R = [S.lambda_mz' S.lambda_s' S.lambda_grand' S.lambda_cousin' S.h_all' S.h_pop' S.h_pathway S.c_R S.h_phantom'];
    R = [R_header' R']';     
    S_latex = latex(R, 2, 3);  % write them also as output to latex file 
    S_latex = mat2cell(S_latex, ones(size(S_latex,1),1));
    savecellfile(S_latex, [schizophrenia_outfile(1:end-4) '_latex_table.txt']); 
end % run schizophrenia data 


ttt = cputime-ttt
