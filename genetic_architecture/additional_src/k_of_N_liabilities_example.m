% A script for computing example of k-of-N liabilities model (Quantitative or binary trait) .
% Give all tables. Set all model parameters
%run_main = 0; % This runs the disease case
compute_power = 0; % compute power to detect loci
plot_power = 0; % plot power figures
compute_QTL=1; % Compute for QTL
plot_QTL = 1; % Plot for QTL
run_binary = 0; % compute deviation for disease model
plot_binary = 0; % plot heritability deviation for binary scale
isoheritability_flag = 0; % 1 means force all QTLs to have the same heritability (heavier). 0 means parameters are the same


AssignGeneralConstants;
ttt = cputime;
K=1; N=3; % # of liabilities and how many are needed to get the disease
freq = 0.05; % minor allele frequency
h_x = 0.999999; % 0.9999999999; % 0.5 % 99999999999; % heritability of EACH liability
num_loci = 20; % num. loci in EACH liability
lambda_mz = 3; % twin risk lambda_mz of disease
%h_x_one_locus =  h_x / num_loci; % 0.025; % assume a given locus contributes this much to the heritability
set_mu_l=0; % we set total disease mu and compute mu_l (otherwise set mu and determine mu_l of a single liability)
if(set_mu_l)
    mu_l = 0.4; % 1/300; % 1/300; % 'prevalence' of one liability of one
    mu = sum(binopdf(K:N, N, mu_l)); % prevalence of disease
else
    mu_vec = [0.001 0.01 0.1]; % This is the 'big official' run for final figure
%    mu_vec = 0.01; % 0.01; % 0.5; % 0.01; % 0.5; % 0.01; % 0.5; % 0.01; % 0.01; % 0.01; % 0.5; % [0.01]; % [0.5]; % [0.01]; % [0.001 0.01 0.1]; %  0.1 0.5]; % Set disease prevalence
    mu_l_vec = zeros(size(mu_vec));
    for i=1:length(mu_vec)
        mu_l_vec(i) = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu_vec(i))), 0, 1); % find mu_l that keeps the prevalence
    end
end
h_shared_env_vec = [0 0.25 0.5 0.75];  % c_R - the fraction of shared environment %  1]; % [0.0 0.5 1]; % How much of the remaining variance is due to shared environment (of twins/sibs etc.)
h_pop_str = {'ACE', 'ADE', 'MZ', 'DZ', 'PO'}; num_h_pop = length(h_pop_str);

test_type_vec = [repmat({'epistasis'}, 1, 10) {'armitage', 'armitage', 'armitage', 'epistasis', 'epistasis', 'epistasis'}];
test_stat_vec = {'probit', 'probit',  'goodness-of-fit',  'goodness-of-fit', ...
    'LRT', 'LRT', 'LRT-analytic', 'gof2d', 'gof2d', 'gof2d-analytic', ...
    'chi-square-analytic', 'chi-square-analytic', 'chi-square', 'probit', 'probit-chi-square', 'probit-analytic'};
sampling_type_vec = [repmat({'full-case-control'}, 1, 10) ...
    repmat({'case-control'}, 1, 6)]; % repmat({'population'}, 1, 6)]; % repmat({'population'}, 1, 6)]; % repmat({'case-control'}, 1, 6)]; %

%    'case-control', 'case-control', 'case-control'}];
test_name_vec = {'pathway_epistasis_probit', 'pathway_epistasis_probit_LT', ...
    'pathway_epistasis_gof', 'pathway_epistasis_gof_LT', ...
    'pathway_epistasis_LRT', 'pathway_epistasis_LRT_LT', 'pathway_epistasis_LRT-analytic',...
    'pathway_epistasis_gof2d', 'pathway_epistasis_gof2d_LT', 'pathway_epistasis_gof2d-analytic', ...
    'one_locus_marginal-analytic', 'one_locus_marginal-analytic', 'one_locus_marginal', ...
    'pairwise_epistasis', 'pairwise_epistasis-chi-square','pairwise_epistasis-analytic'};
power_iters = 500; % power calculation
power_alpha = 5*10^(-8); % genome-wide significance level (marginal test)
power_pairwise_alpha = 0.01 / (0.5*(N*num_loci)*(N*num_loci-1)); % p-value for pairwise test % 5*10^(-6); % more stringent to account for 150^2/2
power_pathway_alpha = 0.01 / (0.5*N*(N-1)); % p-value for pathway test
power_pathway_alpha_vec = [repmat(power_pathway_alpha, 1, 10) ...
    power_alpha power_pathway_alpha power_alpha ...
    power_pairwise_alpha power_pairwise_alpha power_pairwise_alpha];

num_test_stats = length(test_type_vec);
%try_test_inds = [12 13 14 15 16]; % this is single locus and pairwise. 10 is pathway  % [7 10 11 12 13]; % just LRT, gof2d analytic and marginal, and pairwise epistasis
%try_test_inds = [11]; % only single locus (analytic)
try_test_inds = [11 16]; %  16]; % only pairwise epistasis
%%try_test_inds = [10 11 16]; % only analytics. For final figure !!!!

num_mu = length(mu_vec);
power_table_vec = [0.01 0.1 0.25 0.5 0.75 0.9 0.99]; % 'landmark' numbers of power we want to record sample size for
%n_samples_vec = 2*round(logspace(2,7,49)/2); % final, 'official' for paper's figure
n_samples_vec = 2*round(logspace(3,7,999)/2);
%n_samples_vec =  [50000 100000]; % 2*round(logspace(2,5.0,9)/2); % final, 'official' for paper's figure
%[2000:2000:10000];
% 2*round(logspace(2,5,9)/2); %[10^5 3*10^5]; %  2*round(logspace(2,6,21)/2); % 2*round(logspace(2,6,50)/2);
num_n_samples = length(n_samples_vec); % final, 'official' for paper's figure
%h_x_one_locus_vec = logspace(-5,-1.301,49); % final 'official' for paper's figure
h_x_one_locus_vec = 0.0036; % 036; % logspace(-3,-2,19);
%h_x_one_locus_vec = 0.05; %%%% logspace(-5,-1.301,5);
% logspace(-5,-1.301,9); % 0.04; % [0.003 0.01 0.03]; % logspace(-4,-1.301,9); % logspace(-4,-1.301,5); % logspace(-4,-1.301,15); % logspace(-4,-1.301,15); % logspace(-5,-1,45); % heavy run !!!! [1 2 2.5 5 10] ./ 100; % Heritability explained by one locus on liability scale (assuming LT model) %    [0.1 0.5 1 5 10] ./ 100; % [0.005 0.1 0.2 0.5 1 2 5 10 20] ./ 100; % This takes time (linear)  %    logspace(-3,-1,5); % possible heritabilities on X axis
h_x_MLT_one_locus_vec = zeros(num_mu, length(h_x_one_locus_vec));

power_pathway_epistasis_cell = cell(num_test_stats,num_mu);
power_marginal = zeros(num_n_samples, num_mu);

for mu_ind = 1:length(mu_vec)
    mu = mu_vec(mu_ind); mu_l = mu_l_vec(mu_ind);
    % % %     if(run_main) % Run main computation (what example are we using here?)
    % % %         display_parameters = 0;
    % % %         if(display_parameters)% Output numbers we need to fill
    % % %             S
    % % %             mu
    % % %             V = mu*(1-mu)
    % % %             GRR_liability
    % % %             disease_grr_vec
    % % %             p_x_y_is = p_x_y'
    % % %             p_y_given_x_one
    % % %             p_x_z_null = [1-freq freq]' * [1-mu mu]
    % % %             p_x_z
    % % %             p_z_given_x_one
    % % %             p_x_z_diff = p_x_z - p_x_z_null
    % % %             p_x_x_z_disp = [mat2vec(p_x_x_z_zero) mat2vec(p_x_x_z_one)]
    % % %             p_x_x_z_null = [mat2vec(p_x_x_z_zero_null) mat2vec(p_x_x_z_one_null)]
    % % %             p_z_given_x1_x2_one_one = p_z_given_x1_x2(2,2)
    % % %             p_z_given_x1_x2_null_one_one = p_z_given_x1_x2_null(2,2)
    % % %             p_x_x_z_diff = p_x_x_z_disp - p_x_x_z_null
    % % %             h_from_lambda_s = familial_risk_to_heritability(S.lambda_s, 'liability', mu, 0.5)
    % % %             h_from_lambda_mz = familial_risk_to_heritability(S.lambda_MZ, 'liability', mu, 1)
    % % %         end % display parameters
    % % %     end % if run_main
    
    
    if(compute_power) % compute detection power for both main effect and epistasis
        cur_h_x_one_locus_LT = heritability_scale_change_MLT(h_x / num_loci, 1, N, mu, 'LT')/N;
        %    [p_x_y p_x_z p_x_x_z disease_grr_vec] = ...
        %        heritability_to_p_z_x_MLT(cur_h_x_one_locus_LT, freq, K, N, mu);
        
        %    [500:500:4500 5000:5000:200000];
        power_test_type = 'marginal'; % where is this used?
        power_test_stat = 'armitage'; % 'chi-square';
        power_pairwise_test_stat = 'logistic';
        compute_heritability_power_parameters(isoheritability_flag, ...
            power_pathway_alpha_vec, ...
            power_iters, n_samples_vec, ... % power parameters
            h_x, h_x_one_locus_vec, ...
            freq, K, N, mu, mu_ind, mu_vec, num_loci, ... % architecture parameters
            test_name_vec, test_type_vec, test_stat_vec, sampling_type_vec, try_test_inds)        
        
    end % if compute power
end % loop on different prevalences

run_time = cputime - ttt


if(plot_power) % Plot figures showing detection power (marginal and epistasis)
%    star_h_x_vals = [0.36 0.36 0.36]; % frac. variance explained by each SNP for 3x20 loci and ~50% h_pop
        star_h_x_vals = [2.45 0.36 0.36]; % Final! frac. variance explained by each SNP for 3x20 loci and ~50% h_pop
    plot_heritability_power_figures(isoheritability_flag,N, star_h_x_vals);
end % plot power

% main_mu_ind = 3;
if(run_binary) % Compare lambda_s and lambda_mz for disease
    MAX_LAMBDA_MZ = 201; MIN_LAMBDA_MZ = 1.5; % display only values within these ranges
    max_N=10; %10;
    N_vec = 1:max_N; % set independent parameters. Go up to higher N !!!
    
    
    %    mu_vec = 0.01; % [0.001 0.005 0.01 0.05 0.1] % [0.001 0.1]; % try 3 different (shouldn't be very sensitive to mu)
    %    mu_vec = [0.001 0.01]; %  0.05 0.1] % [0.001 0.1]; % try 3 different (shouldn't be very sensitive to mu)
    main_mu_ind = 1; % 3; % which mu to plot as main figure
    % input_lambda_mz_vec = heritability_to_familial_risk(h_x_one_liab, 'liability', mu, 1);
    input_lambda_mz_vec = [2 4]; %  5 10 15]; %  10 15 20]; %  30 40 50 100 200]; % try different values of lambda_mz
    %    input_lambda_mz_vec = [2 20 200]; % try different values of lambda_mz
    
    %    h_x_one_liab_vec = [0.2 0.8]; % [0.1:0.1:0.9]; % need more resolution here
    S_stats = cell(length(mu_vec), length(input_lambda_mz_vec), ...
        length(h_shared_env_vec), max_N, 1);
    lambda_R = S_stats;
    lambda_S_one_liab = S_stats;
    lambda_s_one_liab_vec = zeros(length(mu_vec), length(input_lambda_mz_vec), ...
        length(h_shared_env_vec), max_N, 1);
    lambda_mz_vec = lambda_s_one_liab_vec;
    lambda_s_vec = lambda_s_one_liab_vec;
    disease_mu_vec = lambda_s_one_liab_vec;
    h_liab_loci_vec = lambda_s_one_liab_vec;
    lambda_s_overestimation_vec = lambda_s_one_liab_vec;
    h_phantom_vec = lambda_s_one_liab_vec;
    heritability_overestimation_MZ_vec = lambda_s_one_liab_vec;
    H01 = lambda_s_one_liab_vec;
    h01_loci = lambda_s_one_liab_vec;
    r_MZ01 = lambda_s_one_liab_vec;
    r_DZ01 = lambda_s_one_liab_vec;
    h01_epi = lambda_s_one_liab_vec;
    h01_phantom_vec = lambda_s_one_liab_vec;
    h_liab_pop_vec = cell(num_h_pop,1);
    pi_liab_phantom_vec = cell(num_h_pop,1);
    for i=1:num_h_pop
        h_liab_pop_vec{i} = lambda_s_one_liab_vec;
        pi_liab_phantom_vec{i} = lambda_s_one_liab_vec;
    end
    
        
    
    mu_multidim_vec = lambda_s_one_liab_vec;
    for i=1:length(mu_vec) % 0.04 % [0.001 0.01 0.05 0.1]
        mu = mu_vec(i)
        h_x_one_liab_vec = zeros(1, length(input_lambda_mz_vec));
        if(isoheritability_flag)
            for j=1:length(input_lambda_mz_vec)
                h_x_one_liab_vec(j) = familial_risk_to_heritability(input_lambda_mz_vec(j), 'liability', mu, 1);
            end
        else
            h_x_one_liab_vec = [0.1 0.3 0.5 0.7 0.9]; % big official  run % [0.1 0.3 0.5 0.7 0.9 0.95];
%            h_x_one_liab_vec = [0.5]; % [0.1 0.3 0.5 0.7 0.9]; % [0.1 0.3 0.5 0.7 0.9 0.95];
        end
        
        for j=1:length(h_x_one_liab_vec)
            h_x_one_liab = h_x_one_liab_vec(j) % 0.8 % [0.2 0.5 0.8]
            lambda_mz_one_liab = heritability_to_familial_risk(h_x_one_liab, 'liability', mu, 1);
            
            for i_shared = 1:length(h_shared_env_vec)
                cur_h_shared_env = h_shared_env_vec(i_shared) * (1-h_x_one_liab_vec(j)); % Take Fraction of shared environment
                
                for N=N_vec % loop on different MLT models
                    run_N = N
                    for k=[1]
                        run_k = k
                        sprintf('Run: h_x=%1.2f h_shared=%1.2f mu=%1.2f N=%ld k=%ld', ...
                            h_x_one_liab, cur_h_shared_env, mu, N, k)
                        mu_multidim_vec(i,j,i_shared,N,k) = mu_vec(i);
                        options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance
                        mu_l = fminbnd(@(x) abs(binocdf(k-1, N, x)-(1-mu)), 0, 1, options); % find mu_l that keeps the PREVALENCE
                        switch isoheritability_flag
                            case 0 % just compute disease parameters
                            case 1 % Here fit h_x_one_liab such that total heritability is fixed
                                h_x_one_liab = fminbnd(@(x) abs(lambda_mz_one_liab - compute_k_of_N_liabilities_statistics(...
                                    N, k, mu_l, x, cur_h_shared_env, 1)), 0, 1, options); % compute only for MZ to save time
                        end
                        [lambda_R{i,j,i_shared,N,k} S_stats{i,j,i_shared,N,k}] = ...
                            compute_k_of_N_liabilities_statistics(...
                            N, k, mu_l, h_x_one_liab, cur_h_shared_env, [], h_pop_str); % h_x_one_liab
                        temp_heritability_one_liab = familial_risk_to_heritability(lambda_R{i,j,i_shared,N,k}(1), ...
                            'liability', mu, 1);
                        lambda_S_one_liab{i,j,i_shared,N,k} = ...
                            heritability_to_familial_risk(temp_heritability_one_liab, ...
                            'liability', mu, 0.5);
                        
                        % % %                     % Compute lambda_s over-estimation
                        % % %                     AAA = heritability_to_familial_risk( ...
                        % % %                         familial_risk_to_heritability(A(1), 'liability', mu, 1), ...
                        % % %                         'liability', mu, 0.5)
                        % % %                     (A(2) -AAA) / AAA
                        %                    compute_k_of_N_liabilities_statistics( 1, 1, mu, h_x_one_liab, cur_h_shared_env)
                        
                        lambda_mz_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.lambda_MZ;
                        lambda_s_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.lambda_s;
                        lambda_s_one_liab_vec(i,j,i_shared,N,k) = lambda_S_one_liab{i,j,i_shared,N,k};
                        disease_mu_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.mu;
                        h_liab_loci_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.h_liab_loci;
                        %                        h_liab_pop_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.h_liab_twins;
                        lambda_s_overestimation_vec(i,j,i_shared,N,k) = (lambda_R{i,j,i_shared,N,k}(2) - lambda_S_one_liab{i,j,i_shared,N,k}) / lambda_S_one_liab{i,j,i_shared,N,k};
                        h_phantom_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.h_liab_unexplained_gap;
                        
                        for i_pop = 1:num_h_pop
                            h_liab_pop_vec{i_pop}(i,j,i_shared,N,k) = ...
                                S_stats{i,j,i_shared,N,k}.h_liab_pop(i_pop);
                            pi_liab_phantom_vec{i_pop}(i,j,i_shared,N,k) = ...
                                S_stats{i,j,i_shared,N,k}.pi_liab_phantom(i_pop);
                        end
                        heritability_overestimation_MZ_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.h_liab_unexplained_gap_from_MZ;
                        H01(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.H01;
                        h01_loci(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.h01_loci;
                        r_MZ01(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.r_MZ01;
                        r_DZ01(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.r_DZ01;
                        h01_pop(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.h01_pop;
                        h01_phantom_vec(i,j,i_shared,N,k) = S_stats{i,j,i_shared,N,k}.h01_phantom_vec;
                        
                    end
                end % loop on N
            end % loop on shared environment
        end % loop on heritability
    end % loop on prevalence
    lambda_mz_one_dim_vec = lambda_mz_vec(:);
    mu_one_dim_vec = mu_multidim_vec(:);
    good_inds = intersect(find(lambda_mz_one_dim_vec < MAX_LAMBDA_MZ), ...
        find(lambda_mz_one_dim_vec > MIN_LAMBDA_MZ));
    good_inds = intersect(good_inds, find(mu_one_dim_vec .* lambda_mz_one_dim_vec < 1)); % throw away too large lambda_mz
    
    save(['MLT_data_iso_' num2str(isoheritability_flag) '_debug'], ...
        'S_stats', 'lambda_R', 'h_phantom_vec', ...
        'input_lambda_mz_vec', 'h_liab_loci_vec', 'h_liab_pop_vec', ...
        'lambda_s_vec', 'lambda_s_one_liab_vec', ...
        'lambda_s_overestimation_vec', 'good_inds', 'lambda_mz_vec', 'mu_multidim_vec', ...
        'mu_vec', 'h_x_one_liab_vec', 'disease_mu_vec', 'main_mu_ind', 'max_N', ...
        'h_shared_env_vec', ...
        'MAX_LAMBDA_MZ', 'MIN_LAMBDA_MZ', 'isoheritability_flag', ...
        'H01', 'h01_loci', 'r_DZ01', 'h01_epi', 'h01_phantom_vec', 'pi_liab_phantom_vec');
end % show deviation from LT model for the MLT model (computation)

if(plot_binary)
    load(['MLT_data_iso_' num2str(isoheritability_flag) '_debug']);
    plot_heritability_deviation_figures(['MLT_data_iso_' num2str(isoheritability_flag) '_debug']);
    %    plot_heritability_deviation_figures('MLT_data_iso_1_big');
    %     h_phantom_vec, ...
    %         lambda_s_overestimation_vec, good_inds, lambda_mz_vec, lambda_s_vec, mu_multidim_vec, ...
    %         mu_vec, h_x_one_liab_vec, disease_mu_vec, main_mu_ind, max_N, ...
    %         MAX_LAMBDA_MZ, MIN_LAMBDA_MZ);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Below is QTL Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

operation = 'MIN'; % 'LP'; % 'DIFF'; % 'DIFF'; % 'MAX'; % 'ilogit'; % 'EXP-SUM'; % 'CHI'; % 'MAX'; % 'CHI'; % 'PROD'; % 'MAX';
operation_param = -0.5; % 10.00008;
% operation_param = [3 1];
%operation_param = [-1 0.5];

if(compute_QTL) % Compute statistics for epistasis QTL models
    N_vec = [1 2 5 10 100 1000]; %  100];  % This part isn't heavy - can use many Gaussians
    % % % % % %     if(plot_QTL)
    % % % % % %         plot_max_of_gaussians(N_vec, ['../../common_disease_model/figs/' ...
    % % % % % %             operation '_of_qtls_distribution']);
    % % % % % %     end
    
    main_ind = 1; % 8; % all possible values for h_x
%    h_x_vec = [0.1 0.3 0.5 0.7 0.9]; % big 'official' run
    h_x_vec = [0.3]; 
    num_h = length(h_x_vec);% How much heritability in each liability     %    h_x_vec = [0.1:0.1:0.9]; %  [0.4 0.8];
    
    num_shared_env = length(h_shared_env_vec); % number of different strenghts of shared environments
    max_N = 3; iters = 50000;
    res = 1;
    %    N_vec = res:res:max_N;
    %    N_vec = [N_vec 20 30]; %  50]; % Add more points
    %    N_vec = [1:5]; %  10]; %
 %   N_vec = 1:10; % big official run [1:10]; %  15 20 30 40 50 75 100 200 300 500 1000];
    N_vec = 1:3
    
    num_N = length(N_vec); % which N's to apply
    %     qtl_R = zeros(num_N, 6); % correlation in qtl value for family relatives
    mu = zeros(num_h,num_shared_env,num_N); sigma = zeros(num_h,num_shared_env,num_N);
    h_z = zeros(num_h,num_shared_env,num_N);
    
    h_pop = cell(num_h_pop,1);
    for i=1:num_h_pop % different h_pop esitmations
        h_pop{i} = zeros(num_h,num_shared_env,num_N); % h_pop_ADE=h_pop;
    end
    h_x_output_mat = zeros(num_h,num_shared_env,num_N);
    k_of_N_hist= cell(length(h_x_vec),num_shared_env);
    k_of_N_bins_loc = cell(length(h_x_vec),num_shared_env);
    %    corr_vec = cell(num_N,1); h_z = zeros(num_N,1);
    %    for i=1:num_N % loop on number of gaussians
    N = N_vec(end);
    N_is = N
    compute_mode = 'numeric'; % 'numeric'; % 'simulations'; % 'numeric'; % 'simulations'; % 'simulations'; % 'numeric';
    num_bins = 500;
    for i=1:length(h_x_vec) % number of different heritabilities 
        for j=1:num_shared_env % run for different values of shared environment
            cur_h_shared_env = h_shared_env_vec(j) * (1-h_x_vec(i)); % Take Fraction of shared environment
            [mu(i,j,:) sigma(i,j,:) k_of_N_hist{i,j} k_of_N_bins_loc{i,j} corr_vec{i,j} ...
                h_z(i,j,:) tmp_h_pop h_x_output_mat(i,j,:) qtl_R{i,j} same_inds_prob{i,j}] = ...
                compute_k_of_N_gaussian_statistics( repmat(0, 1, N), repmat(1, 1, N), ...
                h_x_vec(i), cur_h_shared_env, iters, operation, operation_param, ...
                N_vec, num_bins, compute_mode, isoheritability_flag, h_pop_str);
            for i_pop=1:num_h_pop
                h_pop{i_pop}(i,j,:) = tmp_h_pop{i_pop};
            end
        end
    end
    save([operation '_k_of_N_Gaussians_data_iso_' num2str(isoheritability_flag)], ...
        'mu', 'sigma', 'k_of_N_hist', 'k_of_N_bins_loc', 'corr_vec', ...
        'h_z', 'h_pop', 'h_pop_str', 'h_x_vec', 'h_x_output_mat', 'h_shared_env_vec', 'qtl_R', ...
        'same_inds_prob', 'isoheritability_flag', 'main_ind', 'N_vec', 'operation');
    
    % % % % % %     for i=1:num_N
    % % % % % %         [k_of_N_hist{i} k_of_N_bins_loc{i}] = normalize_hist(k_of_N_bins_loc{i}, k_of_N_hist{i}, 1);
    % % % % % %         mean_val = mean_hist(k_of_N_bins_loc{i}, k_of_N_hist{i})
    % % % % % %         std_val = std_hist(k_of_N_bins_loc{i}, k_of_N_hist{i})
    % % % % % %         int_val = integral_hist(k_of_N_bins_loc{i}, k_of_N_hist{i})
    % % % % % %     end
    %        close all
    %end
end % compute QTL parameters

%    h_pop = 2 .* (qtl_R(1,:) - qtl_R(2,:)); % observed in epedimiological study

if(plot_QTL)
    close all; 
    load([operation '_k_of_N_Gaussians_data_iso_' num2str(isoheritability_flag)]);
%    plot_distribution_func_of_gaussians(N_vec, k_of_N_hist, k_of_N_bins_loc); % fig_outfil % just plot the phenotypic distirbution
    close all;
    if(~exist('h_shared_env_vec', 'var'))
        h_shared_env_vec = [0 0.5]; %  1];
        save([operation '_k_of_N_Gaussians_data_iso_' num2str(isoheritability_flag)], ...
            'h_shared_env_vec', '-append');
    end
    %    main_ind = 8; % all possible values for h_x - TEMP!
    if(~exist('h_pop_str', 'var') || isempty(h_pop_str))
        h_pop_str = {'ACE', 'ADE', 'MZ', 'DZ', 'PO'}; num_h_pop = length(h_pop_str);
    end
    clear h_x_shared_env_vec; % Temp! make sure it doesn't exist 
    if(~exist('h_x_shared_env_vec', 'var'))
        h_x_shared_env_vec = h_shared_env_vec; % shared environment within a pathway (relative! c_R!!!) 
        h_shared_env_vec = []; % shared environment in the trait (absolute!!!) 
    end
    plot_heritability_deviation_figures_QTL(h_pop, h_pop_str, h_z, qtl_R, operation, ...
        same_inds_prob, N_vec, main_ind, h_x, h_x_vec, h_x_shared_env_vec, h_shared_env_vec, isoheritability_flag);
end % plot QTL parameters


ttt = cputime - ttt


test_goodness_of_fit=0;
if(test_goodness_of_fit)
    R = rand(n_samples_vec, power_iters); % Test goodness of fit (not chi-square)
    P = 0.0+rand(n_samples_vec, 1).*1;
    %P=0.7;
    X_expected = repmat(P, 1, power_iters);
    X = R < X_expected;
    S = sum((X - X_expected).^2 ./ (X_expected));
    figure; hold on; hist_density(S, 100);
    x_vec = 0:1:n_samples_vec(1);
    %        plot(x_vec, chi2pdf(x_vec, n_samples_vec(1)-1), 'r');
    plot(x_vec, normpdf(x_vec, n_samples_vec*(1-mean(P)), ...
        sqrt (n_samples_vec*mean((1-2.*P).^2.*(1-P)./P))  ), 'm');
    std(S)
    std_should_ve = sqrt (n_samples_vec*mean((1-2.*P).^2.*(1-P)./P))
end % test goodness of fit



% % % % % % % Quick test for heritability. Broad sense should be r_MZ
% % % % % % iters = 100000;
% % % % % % x = randn(iters, 1); % genotype vector
% % % % % % eps = randn(iters, 2); % envirounment
% % % % % %
% % % % % % z1 = (x + 1) .* (eps(:,1)./4 +2) + eps(:,1)./3;
% % % % % % z2 = (x + 1) .* (eps(:,2)./4 +2) + eps(:,2)./3;
% % % % % %
% % % % % % H = corr(x,z1).^2
% % % % % % H_again = corr(x,z2).^2
% % % % % % Twin_Corr = corr(z1, z2).^2
% % % % % %
% % % % % % sqrt(Twin_Corr) - H
% % % % % % sqrt(Twin_Corr) - H_again
% % % % % % (1-Twin_Corr) ./ (1-H)
% % % % % %
% % % % % %
% % % % % %

