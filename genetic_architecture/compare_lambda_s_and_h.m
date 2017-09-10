% Compare sibling risk lambda_s and heritability h: 
% Why can they be different
% function compare_lambda_s_and_h()
% We first do it for a multiplicative model:
% log Pr(z=1|x) = \beta + \sum_i \alpha_i x_i
% L = e^{\alpha_i} so \alpha_i = log(L)

AssignGeneralConstants;
test_simulations = 1; compare_flag = 0; one_locus_model = 0; test_family_to_heritability = 1;
figs_outdir = '../../common_disease_model/figures/';

precision = 3;
L = 2.9; prevalence = 0.025; % 0:0.025:0.025;
RAF = 0.3;
n_vec = [2 4]; n = 2;
architecture_str = 'circuit';

%%%%%%%%%%%%%%%%% First test if simulation and analytical computations agree %%%%%%%%
if(test_simulations)
    for i=1:length(n_vec)
        lambda_s_vec = cell(length(n_vec),1);
        n = n_vec(i);
        f_vec = repmat(RAF, n, 1); GRR_marginal = repmat(L, n, 1);
        for j=1:length(prevalence) % try different prevalence values
            mu = prevalence(j);
            [lambda_s_vec{i} lambda_s(i) lambda_mz(i) h_add(i,j) V_add{i} h(i,j) V(i,j) sib_freq_mat] = ...
                genetic_relative_risk_to_heritability(f_vec, GRR_marginal, mu);
        end
        
        % Compare values of lambda_s obtained by simulations and by analytical computations
        params_struct = [];
        params_struct.z_std = 0; % the environmental st.d.
        params_struct.min_freq = 0;
        params_struct.max_freq = 1;
        params_struct.linear_coef_vec = log(GRR_marginal);
        params_struct.circuit = get_gate(MULTIPLICATIVE, n, params_struct.linear_coef_vec);
        params_struct.circuit_gate_params = cell(n+1,1);
        params_struct.circuit_gate_params{n+1}.linear_coef_vec = params_struct.linear_coef_vec;
        params_struct.circuit_gate_params{n+1}.b = log(prevalence(j)) - ...
            sum(log(f_vec .* (exp(params_struct.linear_coef_vec) - 1) + 1));
        if(~check_architecture_validity(architecture_str, params_struct, n))
            error('Bad Architecture! Not Valid!');
        end
        [family_risk{i} family_tree relative_risk sibs_genotype sibs_freqs] = ...
            compute_architecture_family_risk(architecture_str, ...
            f_vec, params_struct, 500000, ...
            'sampling', 1);
        lambda_s_is = lambda_s(i) % analytically computed lambda_s
        lambda_s_sampling_is = family_risk{i}(end-1) / prevalence(j) % sibling is next-to-last
        
        lambda_mz_is = lambda_mz(i)
        lambda_mz_sampling_is = family_risk{i}(end) / prevalence(j) % twins is last
        
    end
    figure; hold on; % new figure:  family graph plot
    family_size = length(family_tree);
    node_shapes = repmat([0 1 1 0], 1, ceil(family_size/4)); % get female/male shapes (box/circle)
    node_shapes = node_shapes(1:family_size);
    graph_draw(family_tree, ...
        'node_labels', num2str_cell(num2cell(family_risk{i} ./ ...
        prevalence), precision), ...
        'node_shapes', node_shapes); % plot \lambda_R for each family member
    title('Relative risk-odds: person (grey) diseased given family members diseased');
end %%%%%%%%%%%%%%%%% End test if simulation and analytical computations agree %%%%%%%%

if(compare_flag)
    precision = 3;
    L = 1.9; prevalence = [0.001 0.005 0.01 0.05 0.1 0.2];
    RAF = 0.2;
    n_vec = [2 4]; n = 10;
    figure; hold on;
    L_vec = 1:0.1:10;
    f_vec = repmat(RAF, n, 1);
    for j=1:length(prevalence) % try different prevalence values
        mu = prevalence(j);
        if(mu^(1/n) > f_vec(1))
            max_L = (1-f_vec(1)) / (mu^(1/n) - f_vec(1)); % this assumes they're all equal
        else
            max_L = 10000;
        end
        max_i = find(L_vec < max_L, 1, 'last');
        for i=1:max_i % length(L_vec)  % don't take illegal L's
            odds_ratio_marginal = repmat(L_vec(i), n, 1);
            
            params_struct = [];
            params_struct.z_std = 0; % the environmental st.d.
            params_struct.min_freq = 0;
            params_struct.max_freq = 1;
            params_struct.linear_coef_vec = log(GRR_marginal);
            params_struct.circuit = get_gate(MULTIPLICATIVE, n, params_struct.linear_coef_vec);
            params_struct.circuit_gate_params = cell(n+1,1);
            params_struct.circuit_gate_params{n+1}.linear_coef_vec = params_struct.linear_coef_vec;
            params_struct.circuit_gate_params{n+1}.b = log(prevalence(j)) - ...
                sum(log(f_vec .* (exp(params_struct.linear_coef_vec) - 1) + 1));
            if(~check_architecture_validity(architecture_str, params_struct, n))
                error('Bad Architecture! Not Valid!');
            end
            [lambda_s_vec{i} lambda_s(i,j) lambda_mz(i,j) ...
                h_add(i,j) V_add{i} h(i,j) V(i,j) sib_freq_mat] = ...
                genetic_relative_risk_to_heritability(f_vec, GRR_marginal, mu);
        end
        plot(h(1:max_i,j), lambda_mz(1:max_i,j), color_vec(j));
        plot(h(1:max_i,j), lambda_s(1:max_i,j), [color_vec(j) '--']);
    end
    
    
    xlabel('heritability h'); ylabel('familial relative risk \lambda_R');
    title('Heritability vs. family risk for various prevalence values');
    legend_vec = num2str_cell(num2cell(mat_into_vec([prevalence' prevalence']')));
    for j=1:length(legend_vec)
        legend_vec{j} = ['\mu=' legend_vec{j}];
        if(mod(j,2) == 1)
            legend_vec{j} = [legend_vec{j} ' (\lambda_{MZ})'];
        else
            legend_vec{j} = [legend_vec{j} ' (\lambda_{S})'];
        end
    end
    legend(legend_vec, 2);
    my_saveas(gcf, fullfile(figs_outdir, ...
        'heritability_vs_monozigotic_twin_and_sibling_risk_mult_model'), format_fig_vec);
    
    figure; hold on;
    for j=1:length(prevalence)
        plot(lambda_s(:,j), lambda_mz(:,j), color_vec(j));
    end
    plot(lambda_s(:), lambda_s(:).^2, 'r*', 'linewidth', 3);
    xlabel('Sibling risk \lambda_S'); ylabel('Monozigotic twin risk \lambda_{MZ}');
    legend([legend_vec' '\lambda_s^2']', 2);
    my_saveas(gcf, fullfile(figs_outdir, ...
        'sibling_vs_monozigotic_twin_risk_mult_model'), format_fig_vec);
    
    prevalence = [0.001 0.005 0.01 0.05 0.1 0.2];
    figure; hold on; h_vec = 0:0.001:1;
    for j=1:length(prevalence) % Just plot theoretical MZ risk
        plot(h_vec, 1 - h_vec + h_vec ./ prevalence(j), color_vec(j));
    end
    xlabel('heritability h'); ylabel('Monozigotic Twin relative risk \lambda_{MZ}');
    title('Heritability vs. Monozigotic Twin risk for various prevalence values');
    legend_vec = num2str_cell(num2cell(mat_into_vec([prevalence]')));
    for j=1:length(legend_vec)
        legend_vec{j} = ['\mu=' legend_vec{j}];
    end
    legend(legend_vec, 2); %  = num2str_cell(num2cell(mat_into_vec([prevalence]')));
    my_saveas(gcf, fullfile(figs_outdir, 'heritability_vs_monozigotic_twin_risk'), format_fig_vec);
    
end % if compare



if(one_locus_model) % compare all possible effects for one locus for Eliana
    one_locus_model_stats_file = '../../common_disease_model/data/one_locus_model_stats.txt';
    N = 2; f = 0.1; f_vec = [f f]; iters = 1;
    compute_method_flag = 'enumerate'; % 'analytic';
    [x_vec p_x_vec x_ind_vec x_ind_mat] = ...
        initilize_x_vec_constants(N, 0, f_vec);
    
    arch_str_vec = {'dominant', 'recessive', 'additive', 'multiplicative'};
    i=1; z =zeros(2^N, length(arch_str_vec));
    
    relative_sampling_iters = 1000;
    
    for arch_str = arch_str_vec %{'dominant', 'recessive', 'additive', 'multiplicative'}
        eval(['params_struct.gate_type(1) = ' upper(arch_str{1}) ';']);
        params_struct.gate_type(2) = AFFINE;
        params_struct.num_clasues = [1 1];
        params_struct.k_in_clause = [N 1];
        params_struct.gate_params = params_struct.gate_type;
        params_struct.z_std = 0;
        params_struct.levels = 2; % second level does an affine correction (avoid just zeros and ones)
        params_struct.linear_coef_vec = 0.5.*ones(1,N);
        
        switch arch_str{1}
            case 'multiplicative'
                params_struct.b = log(0.1); % for multiplicative architecture
                params_struct.linear_coef_vec(:) = (log(0.9)-params_struct.b)/2;
                params_struct.b(2) = 0; params_struct.a(2) = 1;
            case {'dominant', 'recessive'}
                params_struct.b = 0.1;
                params_struct.b(2) = 0.1; params_struct.a(2) = 0.8; % affine transformation at the end
            otherwise
                params_struct.b = 0; % for additive architecture
                params_struct.b(2) = 0.1; params_struct.a(2) = 0.8; % affine transformation at the end
        end
        
        
        
        [params_struct.circuit params_struct.circuit_gate_params] = ...
            set_architecture_circuit(N, arch_str{1}, params_struct);
        z(:,i) = genetic_architecture(x_vec, 'circuit', params_struct, 1);
        p_x_times_z_vec = p_x_vec .* z(:,i);
        [mu(i) V(i) v_env(i) v_gen(i) mz_twin_risk(i) h(i)] = ...
            compute_architecture_statistics('circuit', ...
            f_vec, params_struct, p_x_vec, p_x_times_z_vec, iters, ...
            compute_method_flag);
        
        [v_marginal GRR_marginal p_z_x_marginal] = ...
            compute_architecture_statistics_marginal('circuit', ...
            f_vec, params_struct, p_x_vec, p_x_times_z_vec, iters, ...
            compute_method_flag, mu(i)); % compute several moments and other stuff for architecture
        h_add2(i) =  sum(V(i) - v_marginal) / V(i);
        [family_risk(:,i) family_tree relative_risk subs_genotype sibs_freq] = ...
            compute_architecture_family_risk('circuit', f_vec', params_struct, ...
            relative_sampling_iters, 'enumerate', 2);
        
        lambda_s(i) = family_risk(end-1) / mu(i);
        
        fit_mode_vec = {'circuit', 'additive', 'additive_naive', ...
            'multiplicative', 'multiplicative_naive', 'logistic', 'liability'}
        for j=1:length(fit_mode_vec)
            fit_mode = fit_mode_vec{j};
            
            fit_mode_arch = str2word('_', fit_mode, 1)
            [alpha beta fitted_params_struct] = fit_architecture_parameters(params_struct, ...
                architecture_str, f_vec, fit_mode, GRR_marginal);
            if(j == 4)
                xxx = 12412412
            end
            z_fitted(:,i) = genetic_architecture(x_vec, fit_mode_arch, fitted_params_struct, 1);
            p_x_times_z_vec_fitted = p_x_vec .* z_fitted(:,i);
            [mu_fitted(i,j) V_fitted(i,j) v_env_fitted(i,j) v_gen_fitted(i,j) ...
                mz_twin_risk_fitted(i,j) h_fitted(i,j)] = ...
                compute_architecture_statistics(fit_mode_arch, ...
                f_vec, fitted_params_struct, p_x_vec, p_x_times_z_vec_fitted, 1, ...
                compute_method_flag);
            fitted_family_risk = ...
                compute_architecture_family_risk(fit_mode_arch, f_vec', fitted_params_struct, ...
                relative_sampling_iters, 'enumerate', 2);
            fitted_lambda = fitted_family_risk ./ mu_fitted(i,j);
            lambda_pa_fitted(i,j) = fitted_lambda(end-2);
            lambda_s_fitted(i,j) = fitted_lambda(end-1);
            lambda_mz_fitted(i,j) = fitted_lambda(end);
            
            if(j == 2) % multiplicative model - another computation
                [lambda_s_vec{i} ...
                    lambda_s_add(i) lambda_mz_add(i) h_add(i) V_add{i} ...
                    lambda_s_mult(i) lambda_mz_mult(i) h_mult(i) V_mult(i) sib_freq_mat] = ...
                    genetic_relative_risk_to_heritability(f_vec', vec2column(GRR_marginal), mu(i)); % exp(alpha)
            end
        end
        lambda_s(i) = family_risk(end-1) / mu(i) % new lambda_s (exact)
        should_be_zero_mz_risk = family_risk(end) / mu(i) - lambda_mz_mult(i)
        i=i+1;
    end % loop on different one-locus 2-allele architectures
    one_locus_labels = {'\rho_{AA}' , '\rho_{AB}', '\rho_{BB}', ...
        '\mu', 'V', 'h', 'h_{add}', 'h_{mult}', ...
        '\lambda_{mz}', '\lambda_{mz}-add', '\lambda_{mz}-mult', ...
        '\lambda_s', '\lambda_s-add', '\lambda_s-mult'};
    h = (V-v_env) ./ V; lambda_mz = heritability_to_mz_twin_risk(h, mu);
    %    lambda_s_add = zeros(length(h),1)-1; % we don't know yet how to compute this !
    
    R = [one_locus_labels'  num2cell([z([1 2 4],:)' mu' V' h' h_add' h_mult' ...
        lambda_mz' lambda_mz_add' lambda_mz_mult' ...
        lambda_s' lambda_s_add' lambda_s_mult'])']';
    R = [['Stat.\Model' arch_str_vec]' R]'
    savecellfile(R, one_locus_model_stats_file, [], 1);
    
    stat_vec = {'h', '\lambda_{mz}', '\lambda_s'};
    one_locus_fit_label = my_tensor_prod(stat_vec, fit_mode_vec); % Try again
    RR = [h_fitted lambda_mz_fitted lambda_s_fitted]';
    RR = [one_locus_fit_label num2cell(RR)]
    
end % if do one locus


if(test_family_to_heritability) % test how people compute heritability
    iters = 1000; arch_str = 'circuit'; n_gates = 4;
    max_generations = 1; f_vec = [0.5 0.5]; % one locus
    circuit = zeros(n_gates);
    circuit(1,3) = 1;  circuit(2,3) = 1; circuit(3,3) = SUM; % % try a linear function AND; % AND; % set 'AND' gate
    circuit(3,4) = 1; circuit(4,4) = AFFINE; % set affine transformation
%    circuit(3,5) = 1; circuit(5,5) = RANDOM; % set random noise
%    circuit(4,6) = 1; circuit(5,6) = 1; circuit(6,6) = SUM; % set 'SUM' gate
    
    params_struct.circuit = circuit; params_struct.linear_coef_vec = [1 1];
    params_struct.circuit_gate_params = cell(n_gates,1); % {'', '', 1, 1};
    params_struct.circuit_gate_params{4}.a = 0.3;
    params_struct.circuit_gate_params{4}.b = 0.2; % affine gate parameters
%     params_struct.circuit_gate_params{5}.a = 0.1;
%     params_struct.circuit_gate_params{5}.b = 0; % random gate parameters
    
    for i=1:2 % simulate two families. In each family the last two are twins
        [family_tree family_genotype_vec] = ...
            simulate_family_genotypes(max_generations, f_vec, iters); % simulate genotype vectors for entire family
        ctr=1;
        for j=[ 1:length(family_tree) length(family_tree)] % simulate phenotype of family members (last one twice to get twins)
            z(ctr+(i-1)*(length(family_tree)+1),:) = ...
                genetic_architecture(family_genotype_vec(:,:,j), arch_str, ...
                params_struct, 1); ctr=ctr+1;
        end
    end
    [h_family h_add_family] = ...
        family_segregation_to_heritability(family_tree, z, 'gaussian') % compute heritability from family
    [h_binary h_add_binary] = ...
        family_segregation_to_heritability(family_tree, z, 'binary') % compute heritability from family
    [h_liability h_add_liability h_add_liability2 thresh mean_liability_given_affected] = ...
        family_segregation_to_heritability(family_tree, z, 'liability') % compute heritability from family
    
    params_struct.z_std = 0;
    [x_vec p_x_vec x_ind_vec x_ind_mat] = ...
        initilize_x_vec_constants(2, 0, f_vec, 'enumerate');
    z_vec = genetic_architecture(x_vec, arch_str, params_struct, 1);
    [mu V v_environment v_genetic mz_twin_risk h_enumerate] = ...
        compute_architecture_statistics(architecture_str, ...
        f_vec, params_struct, p_x_vec', p_x_vec' .* z_vec, 1, 'enumerate')    
    [v_marginal GRR_marginal p_z_x_marginal h_enumerate_add] = ...
        compute_architecture_statistics_marginal(architecture_str, ...
        f_vec, params_struct, p_x_vec', p_x_vec' .* z_vec, 1, 'enumerate', mu) % compute several moments and other stuff for architecture
    
    % Compare two values of h using the conversion formula from Dempster    
    Dempster_h_liability = h_binary * ( 2*pi * mu * (1-mu)) / exp(-thresh^2)  
    diff_liability = Dempster_h_liability - h_liability
    Dempster_h_binary = h_liability * exp(-thresh^2) / ( 2*pi * mu * (1-mu))
    diff_binary = Dempster_h_binary - h_binary 
    
    
    
    [mu V v_env v_gen mz_twin_risk h] = ...
        compute_architecture_statistics(arch_str, ...
        f_vec, params_struct, vec2row(p_x_vec), vec2row(p_x_vec) .* z_vec, 1, ...
        'enumerate');
    [v_marginal GRR_marginal p_z_x_marginal h_add] = ...
        compute_architecture_statistics_marginal(arch_str, ...
        f_vec, params_struct, vec2row(p_x_vec), vec2row(p_x_vec) .* z_vec, 1, ...
        'enumerate', mu);
end
